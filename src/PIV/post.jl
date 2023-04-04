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
a vector of `PranaPass` objects.
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

Base.firstindex(p::PranaData) = firstindex(p.pass)
Base.lastindex(p::PranaData) = lastindex(p.pass)
Base.getindex(p::PranaData, i::Int) = p.pass[i]
Base.length(p::PranaData) = length(p.pass)
Base.iterate(p::PranaData, state=1) = state > length(p) ? nothing : (p[state], state + 1)

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
        hypot(u - um, v - vm) / (r_med + f.ϵ₀) #< f.ϵ
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
