using Statistics, ImageCore, CineFiles

"""
    imadjust(img; qmin=0.01, qmax=0.99)

Adjusts the contrast of an image by scaling the image intensities
to the range [0, 1] using the `qmin` and `qmax` quantiles.
"""
function imadjust(img; qmin=0.01, qmax=0.99)
    vmin, vmax = quantile(vec(img), (qmin, qmax))
    clamp01!(scaleminmax(vmin, vmax).(img))
end

function bgsub(img, bg)
    clamp01!(img .- bg)
end

"""
    tif_load(path)

Loads a TIFF image after temporarily disabling garbage collection,
then re-sets GC to its initial state. This is a workaround for
the 75-90% GC overhead otherwise seen when loading TIFF images.
"""
function tif_load(path)
    initial_gc_state = GC.enable(false)
    img = load(path)
    GC.enable(initial_gc_state)
    img
end

function TSI_bg(bgdata)
    if isfile(bgdata.pathA)
        bgA_img = Gray{Float32}.(tif_load(bgdata.pathA))
    else
        bgA_img = mean(bgdata.bgframes) do i
            tif_load(rawdatadir(bgdata.Date, bgdata.path, TSIname(bgdata.ID, i, "A")))
        end
        save(bgdata.pathA, Gray{N0f16}.(bgA_img))
    end
    if isfile(bgdata.pathB)
        bgB_img = Gray{Float32}.(tif_load(bgdata.pathB))
    else
        bgB_img = mean(bgdata.bgframes) do i
            tif_load(rawdatadir(bgdata.Date, bgdata.path, TSIname(bgdata.ID, i, "B")))
        end
        save(bgdata.pathB, Gray{N0f16}.(bgB_img))
    end
    return bgA_img, bgB_img
end

function phantom_bgsub(runmeta)
    datadict, path = produce_or_load(runmeta, datadir("PLIF", "bgsub"); filename=runname, tag=false) do runmeta
        frames = collect(eval(Meta.parse(runmeta.cine_bgframes)))::Vector{Int}
        path = cinepath(runmeta)
        h = CineFiles.CineHeader(path)
        open(path) do f
            tmp = CineFiles.readframe(f, h, 1)
            bg = mean(frames) do k
                CineFiles.readframe!(f, tmp, h, k)
            end
            ax1, ax2 = axes(bg)
            out = Array{Float32}(undef, reverse(size(bg))..., 200)
            n = first(ax1)+last(ax2)
            for k in 1:200
                CineFiles.readframe!(f, tmp, h, k)
                for i in ax1, j in ax2
                    out[n-j,i,k] = clamp01(tmp[i,j] - bg[i,j])
                end
            end
            Dict("data" => out)
        end
    end
    datadict["data"]
end

function overlapimages(img1::Array{T}, img2::Array{T}) where {T}
    pink = XYZ(Lab(100, 100, -100))
    green = XYZ(Lab(100, -100, 100))
    @. RGB((Float32(img1) * pink) + (Float32(img2) * green))
end

function enhance_TSI(runmeta)
    # Load raw images and background-subtract

    LA_path = rawdatadir(runmeta.Date, runmeta.ID, Vortex.TSI_LA_path(runmeta))
    LB_path = replace(LA_path, "LA" => "LB")
    BGA_path = datadir("PIV", "bg", string(runmeta.TSI_bg_path))
    BGB_path = replace(BGA_path, "LA" => "LB")
    LA, LB = load.((LA_path, LB_path))
    if isfile(BGA_path) && isfile(BGB_path)
        BGA, BGB = load.((BGA_path, BGB_path))

        # Remove high-intensity background pixels from raw images
        BGA_bad = BGA .> 0.3
        BGB_bad = BGB .> 0.3
        LA[BGA_bad] .= 0
        LB[BGB_bad] .= 0
        clamp01!(Gray{Float32}.(LA) .- BGA)
        clamp01!(Gray{Float32}.(LB) .- BGB)
    else
        @warn("No background images found for run $(runname(runmeta)).")
    end

    # Remove maximum-intensity pixels
    LA[LA .> 0.999] .= 0
    LB[LB .> 0.999] .= 0

    # Save intensity-adjusted summary image
    # LA_ceil, LB_ceil = imadjust(LA, qmax=0.999), imadjust(LB, qmax=0.999)
    LA, LB
end

