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
function Base.peek(v::AbstractMatrix)
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
shockruns = filter(r -> !ismissing(r.ptrace_path), runlist)
runmeta = rand(eachrow(shockruns))
pranaraw = PranaData(datadir("PIV", "runs", runname(runmeta)))

u, v, status = vector_replacement(pranaraw, pranaraw.aux["C"][:,:,2]' .< 0.03, (3, 3, 3, 5, 5,))# 7))

GOOD = (status .== PEAK1) .| (status .== PEAK2)
U₀ = hypot.(pranaraw.u, pranaraw.v)
U = hypot.(u, v)
Umin, Umax = quantile(U[GOOD], (0.015, 0.995))
f = Figure()
ax = Axis(f[1:2,1], aspect=DataAspect(), title="Unfiltered vectors")
hm1_hi = heatmap!(ax, U₀; colorrange = (0.2Umax, Umax), lowclip=RGBAf(0,0,0,0), colormap=:inferno)
hm1_lo = heatmap!(ax, U₀, colorrange = (Umin, 0.2Umax), highclip=RGBAf(0,0,0,0))
ax2 = Axis(f[1:2,2], aspect=DataAspect(), title="Spurious vectors replaced with median")
hm2_hi = heatmap!(ax2, U; colorrange = (0.2Umax, Umax), lowclip=RGBAf(0,0,0,0), colormap=:inferno)
hm2_lo = heatmap!(ax2, U, colorrange = (Umin, 0.2Umax), highclip=RGBAf(0,0,0,0))
Colorbar(f[1,3], hm1_hi, label="Post-shock velocity [m/s]")
Colorbar(f[2,3], hm1_lo, label="Pre-shock velocity [m/s]")
linkxaxes!(ax, ax2)
linkyaxes!(ax, ax2)
f

##
f, ax, hm = heatmap(v; colorrange = quantile(v[GOOD], (0.015, 0.995)), axis=(; aspect=DataAspect()))
Colorbar(f[1,2], hm, label="z-velocity [m/s]")
f

##
shockruns = filter(r -> !ismissing(r.ptrace_path), runlist)
runmeta = rand(eachrow(shockruns))
pranaraw = PranaData(datadir("PIV", "runs", runname(runmeta)))

u, v, status = vector_replacement(pranaraw, pranaraw.aux["C"][:,:,2]' .< 0.03, (3, 3, 3, 5, 5,))# 7))
LA, LB = Vortex.enhance_TSI(runmeta)

GOOD = (status .== PEAK1) .| (status .== PEAK2)
##
[imadjust(mapwindow(mean, LB, (33, 33), indices=(8:16:size(LA, 1), (8:16:size(LA, 2))), border="reflect")) .> 0.08; GOOD] |> peek

##
runmeta = rand(eachrow(shockruns))
pranaraw = PranaData(datadir("PIV", "runs", runname(runmeta)))
# LA, LB = Vortex.enhance_TSI(runmeta)
# LIT = imadjust(mapwindow(mean, LB, (33, 33), indices=(8:16:size(LA, 1), (8:16:size(LA, 2))), border="reflect")) .> 0.08
# # peek([LIT; repeatedly(erode ∘ dilate, 3)(LIT)])
# flit = mapwindow(LIT, (9, 9), border="reflect") do w
#     count(w) > 15
# end
BAD = mapwindow(pranaraw.aux["C"][:,:,2]', (9, 9), border="reflect") do w
    median!(w)
end
BAD = mapwindow(BAD, (9, 9), border="reflect") do w
    count(w .< 0.025) > 30
end

u, v, status = vector_replacement(pranaraw, BAD, (3, 3, 3, 5, 5,))# 7))
# @. u[BAD] = 0
# @. v[BAD] = 0

peek([BAD; imadjust(hypot.(u, v))])

##
runmeta = rand(eachrow(runlist))
pranaraw = PranaData(datadir("PIV", "runs", runname(runmeta)))
BAD = mapwindow(pranaraw.aux["C"][:,:,2]', (9, 9), border="reflect") do w
    median!(w)
end
BAD = mapwindow(BAD, (9, 9), border="reflect") do w
    count(w .< 0.025) > 30
end

f = Figure(resolution=(900, 1200))
ax = Axis(f[1, 1], aspect=DataAspect(), height=1100)
sg = SliderGrid(f[2, 1], (label = "Brush size", range=1:30, startvalue=1))
cursorsize = Observable(1)
quality = Menu(f[2, 2], options = zip(["Good", "Okay", "Bad"], 2:-1:0), width=90, tellwidth=true)

on(sg.sliders[1].value) do v
    cursorsize[] = v[]#parse(Int, v[])
end

hm = heatmap!(ax, imadjust(hypot.(pranaraw.u, pranaraw.v)))
# bads = heatmap!(ax, BAD, colormap=[RGBA(0,0,0,0), RGBA(1,0,0,0.5)], colorrange=(0,1))
Vortex.PaintingCanvas(BAD, fig=f, axis=ax, cursorsize=cursorsize, 
    colormap=[RGBA(0,0,0,0), RGBA(1,0,0,0.2)], colorrange=(0,1))
f