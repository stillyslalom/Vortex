using DrWatson
@quickactivate "Vortex"

using Vortex
using CSV, DataFramesMeta
using JLD2
using Dierckx
using LsqFit
using Unitful
using PyThermo
using OrderedCollections
using Roots

##
runlist = loadmeta(m -> !ismissing(m.timings_path) && !ismissing(m.registration_path))
Vortex.add_registrations!(runlist)

runlist.xt = map(eachrow(runlist)) do runmeta
    ismissing(runmeta.ptrace_path) && return missing
    xt = load_xt(runmeta)
end

## Add grid
grids = map(eachrow(runlist)) do runmeta
    phantomscale = runmeta.phantomscale
    xx = range(0, step=phantomscale.mm_per_px, length=576) ./ 1000
    z0 = 0.01425 # distance from end wall of shock tube to tip of vortex generator
    zz = range(stop = phantomscale.origin_mm, step=phantomscale.mm_per_px, length=768) ./ 1000 .- z0
    (x = xx, z = zz)
end
runlist.grid = grids

## Utility functions
core_ix(PLIFcore, lr, state) = first.(values(PLIFcore[lr][state]))
core_iz(PLIFcore, lr, state) = last.(values(PLIFcore[lr][state]))
core_it(PLIFcore, lr, state) = collect(keys(PLIFcore[lr][state]))

## Calculate vortex centroid x0 for one run
runmeta = first(eachrow(runlist))
cores = load(datadir("PLIF", "cores", runname(runmeta)*".jld2"))

function calculate_x0(grid, cores)
    xx, zz = grid
    ismissing(cores) && return missing
    itp_xgrid = Spline1D(eachindex(xx), xx, k=1, bc="extrapolate")
    itp_zgrid = Spline1D(eachindex(zz), zz, k=1, bc="extrapolate")

    @. xcore(z, p) = p[1] + p[2] * z

    x0s = zeros(2)
    for (i, lr) in enumerate(("left", "right"))
        length(core_ix(cores, lr, :preshock)) < 4 && return missing
        x = itp_xgrid(core_ix(cores, lr, :preshock))
        z = itp_zgrid(core_iz(cores, lr, :preshock))
        cf = curve_fit(xcore, z, x, [0.0, 1.0])
        x0s[i] = xcore(0, cf.param)
    end
    x0 = mean(x0s)
end

# f = Figure()
# ax = Axis(f[1, 1], aspect=DataAspect())
# plot!(ax, x, z, label="data")
# lines!(ax, xcore.(z, Ref(cf.param)), z, label="fit")

## Calculate x0 for all runs
runlist.x0 = map(eachrow(runlist)) do runmeta
    corepath = datadir("PLIF", "cores", runname(runmeta)*".jld2")
    cores = !isfile(corepath) ? missing : load(datadir("PLIF", "cores", runname(runmeta)*".jld2"))
    calculate_x0(runmeta.grid, cores)
end

## Find median x-centroid for each registration group
gb = groupby(runlist, :registration_path)
transform!(gb, x -> (; x0_median = median(skipmissing(x.x0))))
unique(runlist.x0_median)

## Shift x0 to median
runlist.grid = map(eachrow(runlist)) do runmeta
    xx, zz = runmeta.grid
    (; xx = xx .- runmeta.x0_median, zz)
end

## Rupture initiation 
# runmeta = eachrow(runlist)[10]

function MST_timing_calc(runmeta)
    sc = MST_state(runmeta)
    Ws = soundspeed(sc.driven)*sc.Ms # shock speed in driven gas
    t_shockprop = 10.61u"inch" / Ws

    timings = read_timings(runmeta)
    rupture_lag = (timings[:Rupture1].Delay + timings[:Rupture1].Shock)*u"μs" + t_shockprop

    t0_Phantom = timings[:Death_Star].Delay - round(Int, ustrip(u"μs", rupture_lag))
    t_Phantom = range(t0_Phantom + 50, length=200, step=50) ./ 1e6

    t_TSI = (timings[:PIV_trig].Delay + runmeta.PIV_delay)*u"μs" - round(u"μs", rupture_lag)
    t_err, i_TSI = findmin(t -> abs(t*u"s" - t_TSI), t_Phantom)
    return (; t₀ = ustrip(u"s", rupture_lag),
              t_Phantom, t_TSI = ustrip(u"s", t_TSI), i_TSI)
end
transform!(runlist, AsTable(:) => ByRow(MST_timing_calc) => [:t₀, :t_Phantom, :t_TSI, :i_TSI])

## Convert cores from indices to physical coordinates
function PLIF_cores(runmeta)
    corepath = datadir("PLIF", "cores", runname(runmeta)*".jld2")
    !isfile(corepath) && return missing
    cores = load(corepath)
    xx, zz = runmeta.grid
    itp_xgrid = Spline1D(eachindex(xx), xx, k=1, bc="extrapolate")
    itp_zgrid = Spline1D(eachindex(zz), zz, k=1, bc="extrapolate")
    corephys = OrderedDict(
        :left => OrderedDict(
            :preshock =>    OrderedDict{Int,NamedTuple{(:t, :x, :z), Tuple{Float64, Float64, Float64}}}(), 
            :postshock =>   OrderedDict{Int,NamedTuple{(:t, :x, :z), Tuple{Float64, Float64, Float64}}}(), 
            :postreshock => OrderedDict{Int,NamedTuple{(:t, :x, :z), Tuple{Float64, Float64, Float64}}}()),
        :right => OrderedDict(
            :preshock =>    OrderedDict{Int,NamedTuple{(:t, :x, :z), Tuple{Float64, Float64, Float64}}}(), 
            :postshock =>   OrderedDict{Int,NamedTuple{(:t, :x, :z), Tuple{Float64, Float64, Float64}}}(), 
            :postreshock => OrderedDict{Int,NamedTuple{(:t, :x, :z), Tuple{Float64, Float64, Float64}}}()))

    for lr in ("left", "right")
        for state in keys(cores[lr])
            ti = core_it(cores, lr, state)
            length(ti) == 0 && continue
            xi = itp_xgrid(core_ix(cores, lr, state))
            zi = itp_zgrid(core_iz(cores, lr, state))
            for (i, ti) in enumerate(ti)
                t = runmeta.t_Phantom[ti]
                x = xi[i]
                z = zi[i]
                corephys[Symbol(lr)][state][ti] = (t=t, x=x, z=z)
            end
        end
    end
    return corephys
end
transform!(runlist, AsTable(:) => ByRow(PLIF_cores) => :cores)

## Shock-vortex interaction z & t (if shocked)
shockruns = filter(r -> !ismissing(r.xt) && !ismissing(r.cine_ID), runlist)
runmeta = first(eachrow(shockruns))

function find_SVI(runmeta)
    missings = (missing, missing, missing, missing)
    ismissing(runmeta.xt) && return missings
    cores = runmeta.cores
    ismissing(cores) && return missings
    corefield(cores, lr, state, field) = getfield.(values(cores[lr][state]), field)

    xx, zz = runmeta.grid

    zfits = Dict{Symbol, Vector{Float64}}() #Vector{Float64}[]
    for lr in (:left, :right)
        t = corefield(cores, lr, :preshock, :t)
        isempty(t) && continue
        x = corefield(cores, lr, :preshock, :x)
        z = corefield(cores, lr, :preshock, :z)

        @. z_c(i, p) = p[1]*(1 - exp(-p[2]*(i - p[3])))
        zfit = curve_fit(z_c, [0; t], [0; z], [0.5, 200, 1e-3])
        # push!(zfits, zfit.param)
        zfits[Symbol(lr)] = zfit.param
    end
    length(zfits) == 0 && return missings

    psfits = Dict{Symbol, Vector{Float64}}() #psfits = Vector{Float64}[]
    for lr in (:left, :right)
        t = corefield(cores, lr, :postshock, :t)
        isempty(t) && continue
        x = corefield(cores, lr, :postshock, :x)
        z = corefield(cores, lr, :postshock, :z)

        @. z_c(i, p) = p[1] + p[2]*i
        zfit = curve_fit(z_c, t, z, Float64[1, -200])
        # push!(psfits, zfit.param)
        psfits[Symbol(lr)] = zfit.param
    end
    length(psfits) == 0 && return missings

    t_SVI = fzero(4e-3) do t
        z_pre = mean(zfit -> zfit[1]*(1 - exp(-zfit[2]*(t - zfit[3]))), values(zfits))
        z_post = mean(zfit -> zfit[1] + zfit[2]*t, values(psfits))
        z_pre - z_post
    end
    z_SVI = mean(zfit -> zfit[1] + zfit[2]*t_SVI, values(psfits))
    return (zfits, psfits, t_SVI, z_SVI)
end

transform!(runlist, AsTable(:) => ByRow(find_SVI) => [:zfits, :psfits, :t_SVI, :z_SVI])

# u₀ = isempty(zfits) ? NaN : mean(zfit -> zfit[1]*zfit[2], values(zfits))
# u = isempty(zfits) ? NaN : mean(zfit -> zfit[1]*zfit[2]*exp(-zfit[2]*t_SVI), values(zfits))
# uₚₛ = isempty(psfits) ? NaN : mean(zfit -> zfit[2], values(psfits))

## Extract and save new values to JLD2
foreach(eachrow(runlist)) do runmeta
    newdata = ntuple2dict((;
        grid = (x = runmeta.grid.xx, z = runmeta.grid.zz),
        t₀ = runmeta.t₀,
        t_Phantom = runmeta.t_Phantom,
        t_TSI = runmeta.t_TSI,
        i_TSI = runmeta.i_TSI,
        cores = runmeta.cores,
        zfits = runmeta.zfits,
        psfits = runmeta.psfits,
        t_SVI = runmeta.t_SVI,
        z_SVI = runmeta.z_SVI,
    ))
    JLD2.save(datadir("master", runname(runmeta)*".jld2"),
        Dict{String, Any}(String(k) => v for (k, v) in newdata))
end
