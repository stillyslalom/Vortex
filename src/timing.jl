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
