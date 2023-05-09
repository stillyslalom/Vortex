using DrWatson
@quickactivate "Vortex"

using Vortex, Glob, DataFramesMeta, Dates, CSV
using PressureTraceXT: PressureTrace, xt, shockspeed, triggerindex, filtertrace!
using GLMakie, Colors
using DSP, ImageFiltering, Statistics, Measurements
using PyThermo, Unitful

## Load run & pressure transducer metadata
runlist = loadmeta()

function PT_meta(path)
    data = DataFrame(:path => glob([r"\d*-\d*-\d*.csv"], path))
    data.Date = map(x -> Dates.Date(splitext(basename(x))[1], "yyyy-mm-dd"), data.path)
    sort!(data, :Date)
    return data
end

channelmap = PT_meta(raw"S:\users\shocktube\hardware_laser_camera_windows_optics_etc\PT_locations\channelmaps")
ptloc = PT_meta(raw"S:\users\shocktube\hardware_laser_camera_windows_optics_etc\PT_locations")

## Define functions to extract and plot XT data
function extract_xt(runmeta, channelmaps, ptlocs)
    ptrace = PressureTrace(rawdatadir(runmeta.Date, runmeta.ID, runmeta.ptrace_path), 
                        last(filter(r -> r.Date < runmeta.Date, channelmaps)).path)
    filtertrace!(ptrace)
    select!(ptrace.data, Not(:PT9))
    xtdata = xt(ptrace,
                last(filter(r -> r.Date < runmeta.Date, ptlocs)).path, 
                :PT3 => 10)
    xtdata.Ws = shockspeed(xtdata, 1)
    (; xtdata, ptrace)
end

function plot_xt(xtdata, ptrace, runmeta)
    ##
    f = Figure()

    ax = Axis(f[1, 1], xlabel="t (s)", ylabel="p (psi)")
    colors = distinguishable_colors(length(names(ptrace.data)), [RGB(0, 0, 0), RGB(1, 1, 1)], dropseed=true)
    for (i, (name, pt)) in enumerate(zip(names(ptrace.data), eachcol(ptrace.data)))
        lines!(ax, ptrace.time, pt, label=name, color=colors[i], linewidth=1)
        vlines!(ax, ptrace.time[triggerindex(ptrace, name => 10)], color=colors[i], linestyle=:dash)
    end
    f[1, 2] = Legend(f[1, 1], ax, valign=:top, halign=:left)

    ax2 = Axis(f[2, 1], xlabel="x (m)", ylabel="t (s)")
    lines!(ax2, Measurements.value.(xtdata.x), xtdata.t, color=RGB(0, 0, 0))
    scatter!(ax2, Measurements.value.(xtdata.x), xtdata.t, color=colors)
    Ws = (xtdata.x[end] - xtdata.x[3])/(xtdata.t[end] - xtdata.t[3])
    text!(ax2, 5, 0; text=
    """
    Wₛ = $(string(Ws))
    Ma = $(string(Ws / 348.4))
    """)

    Label(f[0, :], fontsize = 20, text=string(runname(runmeta), " - ", runmeta.MST_gas, " at ", runmeta.MST_psig, "psig"))
    savefigs && save(plotsdir("xt", runname(runmeta) * "_xt.png"), f)
    f
end

## Extract and plot XT data
savefigs = false
foreach(eachrow(filter(:ptrace_path => !ismissing, runlist))) do r
    try
        xtdata, ptrace = extract_xt(r, channelmap, ptloc)
        savefigs && save(plotsdir("xt", runname(r) * "_xt.png"), plot_xt(xtdata, ptrace, r))
        xtdata.x = Measurements.value.(xtdata.x)
        xtdata.Ws = Measurements.value.(xtdata.Ws)
        CSV.write(datadir("timing", "xt", runname(r) * "_xt.csv"), xtdata)
    catch e
        @warn "Failed to extract xt for $(runname(r))"
    end
end
