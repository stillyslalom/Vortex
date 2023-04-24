using DrWatson
@quickactivate "Vortex"
using Vortex
using Interpolations, BasicInterpolators
using GLMakie
using ImageMorphology
using ImageFiltering
using ImageSegmentation
using StructArrays
using PythonCall

sp_itp = pyimport("scipy.interpolate")

repeatedly(f, n) = ∘(ntuple(_ -> f, n)...)

function Base.peek(v::AbstractMatrix)
    f, ax, hm = heatmap(v; axis=(; aspect=DataAspect()))
    Colorbar(f[1,2], hm)
    f
end

function vorticity(u, v, BAD)
    g₁, g₂ = parent.(Kernel.bickley())
    ImageFiltering.mapwindow(StructArray(; u, v, BAD), (3, 3)) do w
        mapreduce(+, eachindex(w)) do i
            w[i].v * g₁[i] * !w[i].BAD - w[i].u * g₂[i] * !w[i].BAD
        end
    end
end

function manual_vector_validation(runmeta)
    pranaraw = PranaData(datadir("PIV", "runs", runname(runmeta)))
    BAD = mapwindow(pranaraw.aux["C"][:,:,2]', (9, 9), border="reflect") do w
        median!(w)
    end
    BAD = mapwindow(BAD, (9, 9), border="reflect") do w
        count(w .< 0.025) > 30
    end
    
    f = Figure(resolution=(1600, 1200))
    ax = Axis(f[1, 1:2], aspect=DataAspect(), width=700)
    sg = SliderGrid(f[2, 1], (label = "Brush size", range=1:30, startvalue=1))
    cursorsize = Observable(1)
    quality = Menu(f[2, 2], options = zip(["Good", "Okay", "Bad"], 2:-1:0), width=90, tellwidth=true)
    
    on(sg.sliders[1].value) do v
        cursorsize[] = v[]
    end
    
    hm = heatmap!(ax, imadjust(hypot.(pranaraw.u, pranaraw.v)))
    # bads = heatmap!(ax, BAD, colormap=[RGBA(0,0,0,0), RGBA(1,0,0,0.5)], colorrange=(0,1))
    pc = Vortex.PaintingCanvas(BAD, fig=f, axis=ax, cursorsize=cursorsize, 
        colormap=[RGBA(0,0,0,0), RGBA(1,0,0,0.2)], colorrange=(0,1))

    ω = Observable(vorticity(pranaraw.u, pranaraw.v, BAD))
    ωmax = @lift maximum(abs, quantile($(ω)[.!BAD], (0.015, 0.995)))
    ωrange = @lift (-$(ωmax), $(ωmax))
    on(pc.data) do BAD
        ω[] = vorticity(pranaraw.u, pranaraw.v, BAD)
    end

    ax2 = Axis(f[1, 3], aspect=DataAspect(), width=700)
    heatmap!(ax2, ω,
        colorrange = ωrange, colormap=:RdBu)
    linkxaxes!(ax, ax2)
    linkyaxes!(ax, ax2)
    wait(GLMakie.Screen(f.scene)) # wait for plot to close

    Dict("mask"=> .!BAD, "quality"=>quality.selection[])
end

function manual_vector_infill(runmeta)
    pranaraw = PranaData(datadir("PIV", "runs", runname(runmeta)))

    BAD = .!(runmeta.mask)
    u, v, status = vector_infill(pranaraw, BAD, 5)

    f = Figure(resolution=(1600, 1200))
    ax = Axis(f[1, 1:2], aspect=DataAspect(), width=700)
    sg = SliderGrid(f[2, 1], (label = "Brush size", range=1:30, startvalue=1))
    cursorsize = Observable(1)
    qualities = Dict(2 => "Good", 1 => "Okay", 0 => "Bad")
    quality = Menu(f[2, 2], options = zip(["Good", "Okay", "Bad"], 2:-1:0), 
        width=90, tellwidth=true, default=qualities[runmeta.quality])

    on(sg.sliders[1].value) do v
        cursorsize[] = v[]
    end

    uv = Observable(vector_infill(pranaraw, BAD, 5))
    ū = @lift hypot.($(uv).u, $(uv).v)
    ulim = @lift quantile(filter(!isnan, $ū), (0.001, 0.999))

    hm1 = heatmap!(ax, hypot.(pranaraw.u, pranaraw.v),
        colorrange = ulim)

    pc = Vortex.PaintingCanvas(BAD, fig=f, axis=ax, cursorsize=cursorsize, 
        colormap=[RGBA(0,0,0,0), RGBA(1,0,0,0.2)], colorrange=(0,1))

    ax2 = Axis(f[1, 3], aspect=DataAspect(), width=700)
    heatmap!(ax2, ū, colorrange = ulim)
    f[2, 3] = buttongrid = GridLayout(tellwidth=false)
    b = Button(buttongrid[1, 1], label = "Infill")
    on(b.clicks) do _
        uv[] = vector_infill(pranaraw, pc.data[], 5)
    end

    linkxaxes!(ax, ax2)
    linkyaxes!(ax, ax2)
    wait(GLMakie.Screen(f.scene)) # wait for plot to close
    Dict("mask"=> .!BAD, "quality"=>quality.selection[])
end

## Load runs with PIV data
runlist = loadmeta() do meta
    !ismissing(meta.TSI_ID)
end
runlist.runname = runname.(eachrow(runlist))

mask_dicts = map(eachrow(runlist)) do runmeta
    data, path = produce_or_load(manual_vector_validation, runmeta, datadir("PIV", "masks", "initial_masks");
                     filename=runname, tag=false)
    data["runname"] = runname(runmeta)
    data
end

leftjoin!(runlist, DataFrame(mask_dicts), on=:runname)
# filter!(:quality => x -> x != 0, runlist)


## In-fill vectors
mask_dicts = map(eachrow(runlist)) do runmeta
    data, path = produce_or_load(manual_vector_infill, runmeta, datadir("PIV", "masks");
                     filename=runname, tag=false)
    data["runname"] = runname(runmeta)
    data
end
# delete old masks from DataFrame
select!(runlist, Not([:mask, :quality]))
leftjoin!(runlist, DataFrame(mask_dicts), on=:runname)

## Mask finalization =========================================================
# runmeta = rand(eachrow(runlist))
function final_infill(runmeta)
    ni, nj = size(runmeta.mask)
    center_good = component_centroids(Int.(runmeta.mask))[2]
    srg = seeded_region_growing(Gray{N0f8}.(opening(runmeta.mask)), 
        [(CartesianIndex(1, 1), 1),
        (CartesianIndex(round.(Int, center_good)), 2), 
        (CartesianIndex(ni, nj), 1)]).image_indexmap
    unlit = srg .== 1
    uv = vector_infill(PranaData(datadir("PIV", "runs", runname(runmeta))), 
        .!(runmeta.mask), 5; unlit)
    Dict{String, Any}("u"=>uv.u, "v"=>uv.v, "status" => uv.status)
end

map(eachrow(runlist)) do runmeta
    data, path = produce_or_load(final_infill, runmeta, datadir("PIV", "infilled");
                     filename=runname, tag=false)
    data["runname"] = runname(runmeta)
    data
end

##

savefigs = false
foreach(eachrow(runlist)) do runmeta
    ni, nj = size(runmeta.mask)
    center_good = component_centroids(Int.(runmeta.mask))[2]
    srg = seeded_region_growing(Gray{N0f8}.(opening(runmeta.mask)), 
        [(CartesianIndex(1, 1), 1),
        (CartesianIndex(round.(Int, center_good)), 2), 
        (CartesianIndex(ni, nj), 1)]).image_indexmap
    unlit = srg .== 1

    f = Figure(resolution=(1600, 1200))
    ax1 = Axis(f[1, 1:2], aspect=DataAspect(), width=700)
    uv = vector_infill(PranaData(datadir("PIV", "runs", runname(runmeta))), 
        .!(runmeta.mask), 5; unlit)
    heatmap!(ax1, hypot.(uv.u, uv.v))
    ax2 = Axis(f[1, 3], aspect=DataAspect(), width=700)
    heatmap!(ax2, unlit)
    savefigs && save(plotsdir("PIV_masked_magnitude", runname(runmeta)*".png"), f)
end

## Old scratch space ========================================================

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
# BAD = mapwindow(pranaraw.aux["C"][:,:,2]', (9, 9), border="reflect") do w
#     median!(w)
# end
# BAD = mapwindow(BAD, (9, 9), border="reflect") do w
#     count(w .< 0.025) > 30
# end

u, v, status = vector_replacement(pranaraw, BAD, (3, 3, 3, 5, 5,))# 7))
# @. u[BAD] = 0
# @. v[BAD] = 0

peek([BAD; imadjust(hypot.(u, v))])

