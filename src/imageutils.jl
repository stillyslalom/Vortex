using Statistics, ImageCore

function imadjust(img; qmin=0.01, qmax=0.99)
    vmin, vmax = quantile(vec(img), (qmin, qmax))
    clamp01!(scaleminmax(vmin, vmax).(img))
end

function bgsub(img, bg)
    clamp01!(img .- bg)
end

function TSI_bg(bgdata)
    if isfile(bgdata.pathA)
        bgA_img = Gray{Float32}.(load(bgdata.pathA))
    else
        bgA_img = mean(bgframes) do i
            load(rawdatadir(bgdata.Date, bgdata.path, TSIname(bgdata.ID, i, "A")))
        end
        save(bgdata.pathA, Gray{N0f16}.(bgA_img))
    end
    if isfile(bgdata.pathB)
        bgB_img = Gray{Float32}.(load(bgdata.pathB))
    else
        bgB_img = mean(bgframes) do i
            load(rawdatadir(bgdata.Date, bgdata.path, TSIname(bgdata.ID, i, "B")))
        end
        save(bgdata.pathB, Gray{N0f16}.(bgB_img))
    end
    return bgA_img, bgB_img
end

function overlapimages(img1::Array{T}, img2::Array{T}) where {T}
    pink = XYZ(Lab(100, 100, -100))
    green = XYZ(Lab(100, -100, 100))
    @. RGB((Float32(img1) * pink) + (Float32(img2) * green))
end
    