# Postprocessing algorithms for PIV data

## Conversion of Prana data to Julia data structure =================================
struct PranaPass
    x::Vector{Float64}
    y::Vector{Float64}
    u::Matrix{Float64}
    v::Matrix{Float64}
    aux::Dict{String, Any}
end

"""
    PranaPass(matfile::String)

Read a Prana pass file (.mat) and return a `PranaPass` object containing the grid (x & y),
vector field (u & v), and any other data stored in the file (aux). The `aux` field is a
Dict{String,Any} containing all other data stored in the file, including additional correlation
peak displacements (U & V).
"""
function PranaPass(matfile::String)
    vars = matread(matfile)
    X = transpose(pop!(vars, "X"))[:,1]
    Y = transpose(pop!(vars, "Y"))[1,:]
    U = transpose.(eachslice(pop!(vars, "U"), dims=3))
    V = transpose.(eachslice(pop!(vars, "V"), dims=3))
    # TODO: handle case of no extra correlation peak displacements in U & V
    PranaPass(X, Y, collect(U[1]), collect(V[1]), 
        merge!(Dict{String, Any}("U" => U, "V" => V), Dict(string(k) => v for (k, v) in vars)))
end

Base.show(io::IO, p::PranaPass) = print(io, "PranaPass($(length(p.x))×$(length(p.y)))")
Base.size(p::PranaPass) = (length(p.x), length(p.y))
Base.length(p::PranaPass) = length(p.x) * length(p.y)

struct PranaData
    pass::Vector{PranaPass}
end

"""
    PranaData(dir::String)

Read all Prana pass files (.mat) in a directory and return a `PranaData` object containing
a vector of `PranaPass` objects corresponding to each pass. For convenience, the fields of the final pass can
be accessed directly as properties of the `PranaData` object, e.g. `p.x`, `p.y`, `p.u`, `p.v`, `p.aux`.
"""
function PranaData(dir::String)
    pass = PranaPass[]
    for f in readdir(dir)
        if startswith(f, "PIVpass") && endswith(f, ".mat")
            push!(pass, PranaPass(joinpath(dir, f)))
        end
    end
    PranaData(pass)
end

pass(p::PranaData) = getfield(p, :pass)
Base.firstindex(p::PranaData) = firstindex(pass(p))
Base.lastindex(p::PranaData) = lastindex(pass(p))
Base.getindex(p::PranaData, i::Int) = pass(p)[i]
Base.length(p::PranaData) = length(pass(p))
Base.iterate(p::PranaData, state=1) = state > length(p) ? nothing : (p[state], state + 1)
Base.getproperty(p::PranaData, f::Symbol) = getproperty(last(pass(p)), f)

## Spurious vector detection ========================================================

abstract type VectorMedianMeasure end
const VMM = VectorMedianMeasure

# TODO: implement median magnitude
# """
#     MedianMagnitude()

# Select the median vector by magnitude.
# """
# struct MedianMagnitude <: VectorMedianMeasure end
# function (mm::MedianMagnitude)(U)
#     u
# end

struct MedianComponents <: VectorMedianMeasure end
function (mm::MedianComponents)(u, v, δ)
    um = mapwindow(median, u, (δ, δ))
    vm = mapwindow(median, v, (δ, δ))
    return (; u=um, v=vm)
end

abstract type SpuriousVectorDetector end

struct MedianFilter{MM <: VMM} <: SpuriousVectorDetector
    δ::Int
    ϵ::Float64
    mm::MM
    MedianFilter(δ, ϵ, mm::VMM=MedianComponents()) where {VMM <: VectorMedianMeasure} = new{typeof(mm)}(δ, Float64(ϵ), mm)
end

function (f::MedianFilter)(u, v)
    m = f.mm(u, v, f.δ)
    hypot.(u .- m.u, v .- m.v) .< f.ϵ
end


function (f::MedianFilter)(p::PranaPass)
    m = f.mm(p.u, p.v, f.δ)
    map(m.u, m.v, p.u, p.v) do um, vm, u, v
        hypot(u - um, v - vm) < f.ϵ
    end
end

struct NormalizedMedianFilter{MM <: VMM} <: SpuriousVectorDetector
    δ::Int
    ϵ::Float64
    ϵ₀::Float64
    mm::MM
    function NormalizedMedianFilter(δ, ϵ, ϵ₀, mm::VMM=MedianComponents()) where {VMM <: VectorMedianMeasure}
        new{typeof(mm)}(δ, Float64(ϵ), Float64(ϵ₀), mm)
    end
end

function (f::NormalizedMedianFilter)(p::PranaPass)
    m = f.mm(p.u, p.v, f.δ)
    r = map(m.u, m.v, p.u, p.v) do um, vm, u, v
        hypot(u - um, v - vm)
    end
    r_med = mapwindow(median, r, (f.δ, f.δ))
    map(m.u, m.v, p.u, p.v, r_med) do um, vm, u, v, r_med
        hypot(u - um, v - vm) / (r_med + f.ϵ₀) < f.ϵ
    end
end

struct NeighborDifference <: SpuriousVectorDetector
    δ::Int
    ϵ::Float64
end

struct GlobalHistogram <: SpuriousVectorDetector
    N::Int # number of clusters
end

struct DynamicMean <: SpuriousVectorDetector
    δ::Int
    C₁::Float64
    C₂::Float64
end

## Vector replacement ===============================================================
"""
    vector_replacement(u₀, v₀, BAD, r::Int, thresh=10)

Replace spurious vectors in `u₀` and `v₀` with the median of the surrounding `r`×`r` neighborhood.

"""
function vector_replacement(u₀, v₀, BAD, r::Int, thresh=10)
    BAD = BAD .| .!Vortex.MedianFilter(r, thresh)(u₀, v₀)
    u′ = copy(u₀)
    v′ = copy(v₀)
    @. u′[BAD] = NaN
    @. v′[BAD] = NaN
    r₂ = r ÷ 2
    u_clean = mapwindow(u′, (r, r)) do w
        w′ = filter(!isnan, w)
        length(w′) > (r₂) ? median!(w′) : NaN
    end

    v_clean = mapwindow(v′, (r, r)) do w
        w′ = filter(!isnan, w)
        length(w′) > (r₂) ? median!(w′) : NaN
    end

    @. u′[BAD] = u_clean[BAD]
    @. v′[BAD] = v_clean[BAD]
    u′, v′
end

@enum VectorStatus PEAK1 PEAK2 INTERP FAILED=-1

function vector_replacement(pranaraw, BAD, r; thresh=10)
    u₀, v₀ = pranaraw.u, pranaraw.v
    u_init, v_init = copy(u₀), copy(v₀)
    status = Array{VectorStatus}(undef, size(u₀))

    for rᵢ in r
        u₀, v₀ = vector_replacement(u₀, v₀, BAD, rᵢ, thresh)
    end

    u₀[isnan.(u₀)] .= u_init[isnan.(u₀)]
    v₀[isnan.(v₀)] .= v_init[isnan.(v₀)]

    Eval = pranaraw.aux["Eval"][:,:,1]'
    for i in CartesianIndices(u₀)
        if isnan(u₀[i])
            status[i] = FAILED
        elseif u₀[i] != u_init[i]
            status[i] = INTERP
        elseif Eval[i] == 0
            status[i] = PEAK1
        elseif Eval[i] == 1
            status[i] = PEAK2
        else # Eval[i] == 2
            status[i] = FAILED
        end
    end
    u₀, v₀, status
end

const sp_itp = PythonCall.pyimport("scipy.interpolate")

function vector_infill(pranaraw, BAD, r; thresh=50, unlit=falses(size(BAD)))
    BAD .|= .!(Vortex.MedianFilter(r, thresh)(pranaraw.u, pranaraw.v))

    xy_good = collect(reinterpret(reshape, Float64, [Float64.(i.I) for i in findall(.!BAD)])')
    xy_bad = collect(reinterpret(reshape, Float64,  [Float64.(i.I) for i in findall(BAD)])')
    u_spl = sp_itp.LinearNDInterpolator(xy_good, pranaraw.u[.!BAD])
    v_spl = sp_itp.LinearNDInterpolator(xy_good, pranaraw.v[.!BAD])
    u_itp = copy(pranaraw.u)
    v_itp = copy(pranaraw.v)
    u_itp[BAD] .= pyconvert(Vector{Float64}, u_spl(xy_bad))
    v_itp[BAD] .= pyconvert(Vector{Float64}, v_spl(xy_bad))
    u_itp[unlit] .= NaN
    v_itp[unlit] .= NaN

    status = Array{VectorStatus}(undef, size(u_itp))

    Eval = pranaraw.aux["Eval"][:,:,1]'
    for i in CartesianIndices(u_itp)
        if isnan(u_itp[i])
            status[i] = FAILED
        elseif u_itp[i] != pranaraw.u[i]
            status[i] = INTERP
        elseif Eval[i] == 0
            status[i] = PEAK1
        elseif Eval[i] == 1
            status[i] = PEAK2
        else # Eval[i] == 2
            status[i] = FAILED
        end
    end
    (; u=u_itp, v=v_itp, status)
end

function load_infilled_PIV(runmeta)
    d = load(datadir("PIV", "infilled", runname(runmeta)*".jld2"))
    u, v, status = d["u"], d["v"], d["status"]
end