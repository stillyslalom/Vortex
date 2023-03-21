module Vortex

using DrWatson, Reexport
@reexport using DataFramesMeta
@reexport using Statistics
@reexport using ImageCore
using XLSX
using Dates
using CSV
using ImageFiltering, MAT, JLD2

include("pathutils.jl")
include("timing.jl")
include("PIV/post.jl")

export loadmeta, along, select_run

export MedianMagnitude, MedianComponents
export MedianFilter, NormalizedMedianFilter, NeighborDifference, GlobalHistogram, DynamicMean

function loadmeta(f=(_ -> true))
    filter!(f, DataFrame(XLSX.readtable(datadir("meta.xlsx"), "Sheet1")))
end

function along(f::F, A; dims) where {F}
    dropdims(f(A; dims); dims)
end

select_run(runlist, date, ID) = only(filter(r -> r.Date == date && r.ID == ID, runlist))
function select_run(runlist, date_ID)
    date_str, ID = split(date_ID, '_', limit=2)
    select_run(runlist, Date(date_str), ID)
end


end # module Vortex