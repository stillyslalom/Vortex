using DrWatson
@quickactivate "Vortex"
using XLSX, DataFramesMeta, GLMakie

includet(srcdir("imageutils.jl"))
includet(srcdir("pathutils.jl"))
includet(srcdir("timing.jl"))

savefigs = true

## PLIF background subtraction
runlist = DataFrame(XLSX.readtable(datadir("meta.xlsx"), "Sheet1"))

PLIFlist = filter(runlist) do m
    valid = !ismissing(m.cine_ID)
    valid &= !ismissing(m.timings_path)
end

## Write summary plot
foreach(eachrow(PLIFlist)) do PLIFmeta
    plotpath = plotsdir("mean_spanwise_PLIF", runname(PLIFmeta)*".png")
    isfile(plotpath) && return nothing
    timings = read_timings(PLIFmeta)
    # time between trigger (if any) and rupture laser pulse
    rupture_lag = (timings[:Rupture1].Shock + timings[:Rupture1].Delay)
    # First Phantom frame occurs at t0_Phantom after rupture
    t0_Phantom = timings[:Death_Star].Delay - rupture_lag

    cine = phantom_bgsub(PLIFmeta)
    I_avg = dropdims(mean(cine, dims=2), dims=2)
    t_Phantom = range(t0_Phantom, length=size(I_avg, 2), step=50) ./ 1000
    z_Phantom = range(170, step=0.3636, length=size(I_avg, 1))

    f = Figure()
    ax1 = Axis(f[1, 1], title=runname(PLIFmeta), xlabel="t [ms]", ylabel="z [mm]")
    hm = heatmap!(ax1, t_Phantom, z_Phantom, rotr90(I_avg))
    savefigs && save(plotpath, f)
    f
end

##