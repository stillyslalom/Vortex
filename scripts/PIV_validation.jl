using DrWatson
@quickactivate "Vortex"
using Vortex
using Interpolations, BasicInterpolators
using GLMakie
using ImageMorphology
using ImageFiltering
using StructArrays


repeatedly(f, n) = ∘(ntuple(_ -> f, n)...)

# Load runs with PIV data
runlist = loadmeta() do meta
    !ismissing(meta.TSI_ID)
end

## Select & load a random run
runmeta = rand(eachrow(runlist))
pranaraw = PranaData(datadir("PIV", "runs", runname(runmeta)))
function peek(v)
    f, ax, hm = heatmap(v; axis=(; aspect=DataAspect()))
    Colorbar(f[1,2], hm)
    f
end

BAD = pranaraw[end].aux["C"][:,:,2]' .< 0.025
BAD .&= .!Vortex.MedianFilter(5, 10)(pranaraw[3])
peek(BAD)
# (repeatedly(erode ∘ dilate, 10))(BAD) |> peek
# peek(pranaraw[end].aux["V"][2])
# peek(MedianFilter(3, 10)(pranaraw[end]))
u′ = copy(pranaraw[end].u)
@. u′[BAD] = NaN
@time u_clean = mapwindow(u′, (3,3)) do w
    w′ = filter(!isnan, w)
    length(w′) > 5 ? median!(w′) : NaN
end
 
@. u′[BAD] = u_clean[BAD]
peek([pranaraw.u; u′])
##

function vector_replacement(u₀, v₀, BAD, r::Int)
    BAD = BAD .& .!Vortex.MedianFilter(r, 10)(u₀, v₀)
    u′ = copy(u₀)
    v′ = copy(v₀)
    @. u′[BAD] = NaN
    @. v′[BAD] = NaN
    r₂ = r ÷ 2
    u_clean = mapwindow(u′, (r, r)) do w
        w′ = filter(!isnan, w)
        length(w′) > (r₂) ? median!(w′) : NaN
    end
    # u_clean = imfilter(u_clean, KernelFactors.gaussian((1, 1), (3, 3)))

    v_clean = mapwindow(v′, (r, r)) do w
        w′ = filter(!isnan, w)
        length(w′) > (r₂) ? median!(w′) : NaN
    end
    # v_clean = imfilter(v_clean, KernelFactors.gaussian((1, 1), (3, 3)))

    @. u′[BAD] = u_clean[BAD]
    @. v′[BAD] = v_clean[BAD]
    u′, v′
end

function vector_replacement(u₀, v₀, BAD, r)
    u_init, v_init = copy(u₀), copy(v₀)
    for rᵢ in r
        u₀, v₀ = vector_replacement(u₀, v₀, BAD, rᵢ)
    end
    u₀[isnan.(u₀)] .= u_init[isnan.(u₀)]
    v₀[isnan.(v₀)] .= v_init[isnan.(v₀)]
    u₀, v₀
end

##
runmeta = rand(eachrow(runlist))
pranaraw = PranaData(datadir("PIV", "runs", runname(runmeta)))

u, v = vector_replacement(pranaraw.u, pranaraw.v, pranaraw.aux["C"][:,:,2]' .< 0.03, (3, 3, 3, 5, 5, 7))
peek([pranaraw.u; u])