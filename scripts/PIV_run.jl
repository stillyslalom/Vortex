using DrWatson
@quickactivate "Vortex"

using XLSX, DataFramesMeta, TOML
using GLMakie, UnPack, ImageFiltering

# Here you may include files from the source directory
includet(srcdir("PIV", "prana.jl"))
includet(srcdir("imageutils.jl"))
includet(srcdir("pathutils.jl"))

## Load run metadata
runlist = DataFrame(XLSX.readtable(datadir("meta.xlsx"), "Sheet1"))
runlist.outdir = [datadir("PIV", "runs", runname(runmeta)) for runmeta in eachrow(runlist)]
runlist.TSI_LA_path = [TSIname(runmeta.TSI_ID, runmeta.TSI_idx, 'A') for runmeta in eachrow(runlist)]
runlist.validpaths = [validatepaths(runmeta) for runmeta in eachrow(runlist)]
runmeta = last(runlist)

## Set up PIV runs
PIVlist = filter(runlist) do m
    # if PIV.toml file exists in output directory, only run if m.PIV_rerun is true.
    runjob = ispath(joinpath(m.outdir, "PIV.toml")) ? m.PIV_rerun : true
    runjob &= !ismissing(m.TSI_idx) # only run if TSI data isn't missing
    runjob &= m.validpaths # only run if all paths are valid
    return runjob
end

## Launch jobs asychronously on a worker process

foreach(eachrow(PIVlist)) do runmeta
    if runmeta.PIV_rerun && ("PIV_rerun.toml" in readdir(runmeta.outdir))
        cfg_m = setupprana(runmeta, toml2config(joinpath(runmeta.outdir, "PIV_rerun.toml")))
    else
        cfg_m = setupprana(runmeta)
    end

    try
        runprana(runmeta, cfg_m)
    catch 
        e
    end
end

##

foreach(eachrow(runlist)) do m
    runjob = !ismissing(m.TSI_idx) # only run if TSI data isn't missing
    runjob &= m.validpaths # only run if all paths are valid
    runjob && pranasummaryplot(m)
end
