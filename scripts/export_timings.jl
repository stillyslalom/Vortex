using DrWatson
@quickactivate "Vortex"
using Vortex

##
runlist = loadmeta() do m
    !ismissing(m.timings_path)
end

##
for m in eachrow(runlist)
    srcpath = rawdatadir(m.Date, m.ID, m.timings_path)
    dstpath = datadir("timing", "tsv", join((string(m.Date), m.ID, m.timings_path), '_'))
    isfile(srcpath) && !isfile(dstpath) && cp(srcpath, dstpath)
end

## test
map(read_timings, eachrow(runlist))