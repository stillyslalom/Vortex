using DrWatson
@quickactivate "Vortex"
using XLSX, DataFramesMeta

includet(srcdir("registration.jl"))
includet(srcdir("imageutils.jl"))
includet(srcdir("pathutils.jl"))

##
targetlist = DataFrame(XLSX.readtable(datadir("targets.xlsx"), "Sheet1"))

targetdata = first(targetlist)

targets = load_targets(targetdata)
run_cpselect(targets)
