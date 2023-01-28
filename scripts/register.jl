using DrWatson
@quickactivate "Vortex"
using XLSX, DataFramesMeta, GLMakie
using ImageTransformations, CoordinateTransformations, StaticArrays

includet(srcdir("registration.jl"))
includet(srcdir("imageutils.jl"))

##
targetlist = DataFrame(XLSX.readtable(datadir("targets.xlsx"), "Sheet1"))

targetdata = first(targetlist)

targets = load_targets(targetdata)
tf = produce_or_load(targetdata, datadir("calibrations");
                     filename=runname, suffix="toml", tag=false) do targetdata
    targets = load_targets(targetdata)
    tf = run_cpselect(targets, mkdir(datadir("calibrations", runname(targetdata))))


run_cpselect(targets)
T = copy(tf)
##
M = LinearMap(SMatrix{3,3,Float32}(T'))
push1(x) = CoordinateTransformations.push(x, 1)

tform = PerspectiveMap() ∘ M ∘ push1
itform = PerspectiveMap() ∘ inv(M) ∘ push1

# TSI_warp = warp(padarray(targets.TSI, Fill(1, (30, 30))), itform, axes(targets.Phantom)) |> imadjust
TSI_warp = warp(targets.TSI, itform, axes(targets.Phantom)) |> imadjust
tgts = overlapimages(Gray{Float64}.(targets.Phantom), TSI_warp)

f, ax, img = image(collect(tgts'); axis=(aspect=DataAspect(), yreversed=true))
# list corners of TSI
corners = [SVector(1, 1), SVector(size(targets.TSI, 1), 1), SVector(size(targets.TSI, 1), size(targets.TSI, 2)), SVector(1, size(targets.TSI, 2))]
push!(corners, first(corners))
ij_to_xy(c) = reverse(tform(c))
horiz_mid = [mean((corners[1], corners[2])), mean((corners[3], corners[4]))]
vert_mid = [mean((corners[1], corners[4])), mean((corners[2], corners[3]))]
scatter!(ax,ij_to_xy.(corners), color=:red, linewidth=2)
lines!(ax, ij_to_xy.(corners), color=:red, linewidth=2)
lines!(ax, ij_to_xy.(horiz_mid), color=:red, linewidth=1)
lines!(ax, ij_to_xy.(vert_mid), color=:red, linewidth=1)

f
##
p = phantomscale(targets.Phantom)
##
px_per_inch = abs(dist2d(p.point1.coords, p.point2.coords) / 
                        (p.point1.pos - p.point2.pos))*u"inch^-1"
mm_per_px = ustrip(u"mm", inv(px_per_inch))
midpoint_px = mean((p.point1.coords, p.point2.coords))
midpoint_mm = ustrip(u"mm", mean((p.point1.pos, p.point2.pos)) * u"inch")
origin_mm = midpoint_mm + midpoint_px[2]*mm_per_px
phantom_z = range(origin_mm, length=size(targets.Phantom, 1), step=-mm_per_px)
phantom_x = range(0, length=size(targets.Phantom, 2), step=mm_per_px)

image(phantom_x, phantom_z, rotr90(targets.Phantom); axis=(aspect=DataAspect(),))
##
@time Phantom_warp = warp(targets.Phantom, tform, axes(targets.TSI)) |> imadjust
tgts = overlapimages(Gray{Float64}.(targets.TSI) |> imadjust, Phantom_warp)




