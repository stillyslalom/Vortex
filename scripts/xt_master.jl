using DrWatson
@quickactivate "Vortex"

using Vortex
using CSV, DataFramesMeta
using JLD2
using Dierckx
using LsqFit
using Unitful
using PyThermo

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

# t_Phantom = range(t0_Phantom + 50, length=size(cine, 3), step=50) ./ 1e6

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
runmeta = eachrow(runlist)[40]
sc = MST_state(runmeta)
Ws = soundspeed(sc.driven)*sc.Ms # shock speed in driven gas
t_shockprop = 10.61u"inch" / Ws

timings = read_timings(runmeta)
rupture_lag = (timings[:Rupture1].Delay + timings[:Rupture1].Shock)*u"μs" + t_shockprop

t0_Phantom = timings[:Death_Star].Delay - round(Int, ustrip(u"μs", rupture_lag))
t_Phantom = range(t0_Phantom + 50, length=200, step=50) ./ 1e6

t_TSI = (timings[:PIV_trig].Delay + timings[:PIV_trig].Shock + runmeta.PIV_delay)*u"μs" - round(u"μs", rupture_lag)
findmin(t -> abs(t*u"s" - t_TSI), t_Phantom)