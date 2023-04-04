using DrWatson
@quickactivate "Vortex"
using Vortex
using Interpolations, BasicInterpolators
using GLMakie
using ImageMorphology
using ImageFiltering

repeatedly(f, n) = ∘(ntuple(_ -> f, n)...)

# Load runs with PIV data
runlist = loadmeta() do meta
    !ismissing(meta.TSI_ID)
end

## Select & load a random run
# runmeta = rand(eachrow(runlist))
pranaraw = PranaData(datadir("PIV", "runs", runname(runmeta)))
function peek(v)
    f, ax, hm = heatmap(v; axis=(; aspect=DataAspect()))
    Colorbar(f[1,2], hm)
    f
end

BAD = pranaraw[end].aux["C"][:,:,1]' .< 0.025
(repeatedly(erode ∘ dilate, 10))(BAD) |> peek
# peek(pranaraw[end].aux["V"][2])
# peek(MedianFilter(3, 10)(pranaraw[end]))