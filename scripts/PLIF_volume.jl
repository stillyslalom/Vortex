using DrWatson
@quickactivate "Vortex"
using Vortex
using CairoMakie
using JLD2
using LaTeXStrings
using ImageFiltering
using Statistics

update_theme!(Theme(fonts = (; regular = "New Roman", bold = "New Roman Bold")))
##

runlist = loadmeta() do m
    all(!ismissing, (m.TSI_ID, m.registration_path, m.cine_ID, m.timings_path))
end

masters = map(eachrow(runlist)) do runmeta
    data = JLD2.load(datadir("master", runname(runmeta)*".jld2"))
    data["runname"] = runname(runmeta)
    data
end |> DataFrame

runlist.runname = runname.(eachrow(runlist))
leftjoin!(runlist, masters, on = :runname)

##

meta = first(eachrow(runlist))
PLIF = phantom_bgsub(meta)

##
frame = PLIF[:,:,19]

# I_bg = median!(filter(!=(0), frame))
# I_thresh = mean(quantile(filter(>(2I_bg), frame), (0.1, 0.90000000000000000000000000000000000000000000000000000000000)))
# (frame .> I_thresh) |> heatmap


mframe = mapwindow(median, frame, (11, 11))

(mframe .> mean(quantile(vec(mframe), (0.05, 0.95))) .* 1.0) |> heatmap
