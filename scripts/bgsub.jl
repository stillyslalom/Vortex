using DrWatson
@quickactivate "Vortex"
using XLSX, DataFramesMeta

includet(srcdir("imageutils.jl"))
includet(srcdir("pathutils.jl"))

## Background subtraction
bglist = DataFrame(XLSX.readtable(datadir("backgrounds.xlsx"), "Sheet1"))
@rtransform! bglist @astable begin
    :bgframes = collect(eval(Meta.parse(:frames)))::Vector{Int}
    :basename = join((:Date, :path, :ID, first(:bgframes), last(:bgframes)), '-')
    :pathA = datadir("PIV", "bg", :basename*"LA.tif")
    :pathB = datadir("PIV", "bg", :basename*"LB.tif")
end
bgdata = first(bglist)

TSI_bg(bgdata)
