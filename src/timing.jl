export read_timings, load_xt

function read_timings(meta)
    timings = CSV.read(rawdatadir(meta.Date, meta.ID, meta.timings_path), DataFrame, 
        normalizenames=true)
    timings.Channel = @. Symbol(replace(strip(timings.Channel), isspace => '_'))
    @chain timings begin
        groupby("Channel")
        filter(g -> size(g, 1) == 1, _)
        [only(r) for r in _]
        Dict(r.Channel => NamedTuple(r[Not(:Channel)]) for r in _)
    end
end

function load_xt(runmeta)
    CSV.read(datadir("timing", "xt", runname(runmeta) * "_xt.csv"), DataFrame)
end

"""
    phantom_timing(runmeta, timings, cine_nframes)

Find the index of the cine frame closest to the PIV imaging time, as well as
the discrepancy between the two times.
"""
function phantom_timing(runmeta, timings, cine_nframes)
    t_TSI = (timings[:PIV_trig].Delay + runmeta.PIV_delay) / 1e6
    t0_Phantom = timings[:Death_Star].Delay + 50
    t_Phantom = range(t0_Phantom, length=cine_nframes, step=50) / 1e6
    t_err, cine_PIV_idx = findmin(t -> abs(t - t_TSI), t_Phantom)
    return cine_PIV_idx, t_err
end

"""
    MST_state(MST_gas, MST_psig)

Calculate the vortex ring generator gas states for a given run.
"""
function MST_state(MST_gas, MST_psig)
    p_ace = Species("acetone", T = 22u"°C").Psat*u"Pa"
    p_total = (MST_psig + 14.5)*u"psi" |> u"Pa"
    χ_ace = p_ace/p_total
    MSTgas = Mixture([MST_gas => 1-χ_ace, "acetone" => χ_ace])
    Ma = fzero(1.31) do M
        pressure(PyThermo.ShockTube.driverpressure(MSTgas, MSTgas, M)) - p_total
    end
    sc = PyThermo.ShockTube.shockcalc(MSTgas, MSTgas, Ma)
end

"""
    MST_state(runmeta)

Calculate the vortex ring generator gas states for a given run.
"""
MST_state(runmeta) = MST_state(runmeta.MST_gas, runmeta.MST_psig)
