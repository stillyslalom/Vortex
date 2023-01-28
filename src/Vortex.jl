module Vortex

using DrWatson, XLSX, DataFramesMeta

include("pathutils.jl")
include("timing.jl")

export loadmeta

function loadmeta()
    DataFrame(XLSX.readtable(datadir("meta.xlsx"), "Sheet1"))
end

end # module Vortex