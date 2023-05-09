module Vortex

using DrWatson, Reexport
@reexport using DataFramesMeta
@reexport using Statistics
@reexport using ImageCore
using XLSX
using Dates
using CSV
using ImageFiltering, MAT, JLD2
using ImageTransformations, StaticArrays, CoordinateTransformations, Interpolations
using StructArrays
using CineFiles
using PythonCall
using GLMakie
using PyThermo
using Unitful
using Roots: fzero

update_theme!(Theme(fonts = (; regular = "New Roman", bold = "New Roman Bold")))

include("pathutils.jl")
include("imageutils.jl")
include("timing.jl")
include("PIV/post.jl")
include("registration.jl")
include("painting.jl")
include("fluiddyn.jl")

export loadmeta, along, select_run
export imadjust, phantom_bgsub, overlapimages

export MedianMagnitude, MedianComponents
export MedianFilter, NormalizedMedianFilter, NeighborDifference, GlobalHistogram, DynamicMean
export PranaData, PranaPass, vector_replacement, VectorStatus, vector_infill
export PEAK1, PEAK2, INTERP, FAILED # vector status enum
export push1
export MST_state
export vorticity

function loadmeta(f=(_ -> true))
    filter!(f, DataFrame(XLSX.readtable(datadir("meta.xlsx"), "Sheet1")))
end

function along(f::F, A; dims) where {F}
    dropdims(f(A; dims); dims)
end

select_run(runlist, date, ID) = only(filter(r -> r.Date == date && r.ID == ID, runlist))

"""
    select_run(runlist, date_ID::String)

Select a run from `runlist` by its date and ID, given as a string in the format
"YYYY-mm-dd_ID".
"""
function select_run(runlist, date_ID::String)
    date_str, ID = split(date_ID, '_', limit=2)
    select_run(runlist, Date(date_str), ID)
end


end # module Vortex