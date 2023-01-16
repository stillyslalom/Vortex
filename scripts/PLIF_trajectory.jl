using DrWatson
@quickactivate "Vortex"
using XLSX, DataFramesMeta, JLD2, GLMakie, LsqFit, ImageFiltering, Peaks, CineFiles

includet(srcdir("imageutils.jl"))
includet(srcdir("pathutils.jl"))

savefigs = true

## PLIF background subtraction
runlist = DataFrame(XLSX.readtable(datadir("meta.xlsx"), "Sheet1"))

PLIFlist = filter(runlist) do m
    valid = !ismissing(m.cine_ID)
end

PLIFmeta = first(PLIFlist)

##
# begin; ydata = eachcol(I_avg)[10]
# I_max = map(eachcol(I_avg)) do ydata
#     xdata = eachindex(ydata)

#     yfilt = imfilter(reverse(ydata), KernelFactors.gaussian(10))
#     pks, vals = findmaxima(yfilt, 10)
#     pks, proms = peakproms!(pks, yfilt)
#     pks, ws, ls, rs = peakwidths!(pks, yfilt, proms; relheight=0.5)
#     peakvol = map(zip(pks, vals, ls, rs)) do (pk, v, l, r)
#         v*(r - l)
#     end
#     # return largest two peaks by volume
#     imax = partialsortperm(peakvol, 1:min(2, length(peakvol)), rev=true)
#     length(peakvol) == 0 && return (0, 0)
#     length(peakvol) == 1 && return fill(pks[imax[1]], 2)
#     return pks[imax]
# end


# I_max1, I_max2 = first.(I_max), last.(I_max)
# I_max_clean = map(zip(mapwindow(median, first.(I_max), (3,)), I_max1, I_max2)) do (xmed, x1, x2)
#     abs(x1 - xmed) < abs(x2 - xmed) ? x1 : x2
# end

# ##
# let
#     f = Figure()
#     ax1 = Axis(f[1, 1])
#     sl = Slider(f[2, :], range=1:200)
#     I_avg = reduce(hcat, [Float32.(mean(rotl90(cf[i] - bg), dims=2)) for i = 1:200])
#     hm = heatmap!(ax1, rotr90(I_avg))
#     idx = sl.value
#     vlines!(ax1, @lift([$idx]))
#     scatter!(ax1, I_max_clean, marker=:cross, color=:red, markersize=10)
#     xlims!(ax1, 1, 60)
#     ax2 = Axis(f[1, 2])
#     xlims!(ax2, (0, maximum(maximum, I_avg)))
#     lines!(ax2, @lift(reverse(vec(I_avg[:, $idx] - median(I_avg, dims=2)))), eachindex(I_avg[:, 1]))
#     lines!(ax2, reverse(vec(median(I_avg, dims=2))), axes(I_avg, 1))
#     f
# end

