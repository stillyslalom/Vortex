using DrWatson
@quickactivate "Vortex"
using Vortex
using CairoMakie
using JLD2
using LaTeXStrings
using Printf

update_theme!(Theme(fonts = (; regular = "New Roman", bold = "New Roman Bold")))
##

runlist = loadmeta() do m
    !ismissing(m.TSI_ID) && !ismissing(m.registration_path) && !ismissing(m.cine_ID) && !ismissing(m.timings_path)
end

masters = map(eachrow(runlist)) do runmeta
    data = JLD2.load(datadir("master", runname(runmeta)*".jld2"))
    data["runname"] = runname(runmeta)
    data
end |> DataFrame

runlist.runname = runname.(eachrow(runlist))
leftjoin!(runlist, masters, on = :runname)

runlist.PIV = [load(datadir("PIV", "registered", r.runname*".jld2")) for r in eachrow(runlist)]
runlist.PIVcore = [load(datadir("PIV", "cores", r.runname*".jld2")) for r in eachrow(runlist)]
runlist.PIV_quality = map(eachrow(runlist)) do runmeta
    maskpath = datadir("PIV", "masks", runname(runmeta)*".jld2")
    !isfile(maskpath) && return missing
    load(maskpath, "quality")
end
##
