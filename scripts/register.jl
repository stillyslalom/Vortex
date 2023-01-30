using DrWatson
@quickactivate "Vortex"
using Vortex
using XLSX, DataFramesMeta, GLMakie
using ImageTransformations, CoordinateTransformations, StaticArrays

includet(srcdir("registration.jl"))
includet(srcdir("imageutils.jl"))

##
targetlist = DataFrame(XLSX.readtable(datadir("targets.xlsx"), "Sheet1"))

targetdata = eachrow(targetlist)[4]

calibration_dicts = map(eachrow(targetlist)) do targetdata
    data, path = produce_or_load(targetdata, datadir("calibrations");
                     filename=runname, tag=false) do targetdata
        targets = load_targets(targetdata)
        calibration_plotdir = datadir("calibrations", runname(targetdata))
        # remove & re-make calibration plot directory if it already exists
        isdir(calibration_plotdir)  && rm(calibration_plotdir; recursive=true)
        mkdir(calibration_plotdir) 

        tf = run_cpselect(targets, calibration_plotdir)

        phantomscalepoints = selectphantomscale(targets.Phantom)
        phantomscale = calcphantomscale(phantomscalepoints...)
        Dict("transform"=>tf, "phantomscale"=>phantomscale)
    end
    data["path"] = path
    data
end

calibrations = hcat(targetlist, DataFrame(calibration_dicts))
##
foreach(eachrow(calibrations)) do calib
    targets = load_targets(calib, datadir("calibrations", runname(calib)))
    M = LinearMap(SMatrix{3,3,Float32}(calib.transform'))
    push1(x) = CoordinateTransformations.push(x, 1)

    tform = PerspectiveMap() ∘ M ∘ push1
    itform = PerspectiveMap() ∘ inv(M) ∘ push1

    TSI_warp = warp(targets.TSI, itform, axes(targets.Phantom))
    tgts = overlapimages(targets.Phantom, TSI_warp)

    f = Figure(resolution=(950, 1200))
    ax = Axis(f[1, 1], aspect=DataAspect(), yreversed=true, 
            xlabel="jᵗʰ Phantom pixel", ylabel="iᵗʰ Phantom pixel",
            title = "Phantom and TSI images (warped to Phantom coordinates) for $(runname(calib))")
    img = image!(ax, collect(tgts'))
    # list corners of TSI
    corners = [SVector(1, 1), SVector(size(targets.TSI, 1), 1), SVector(size(targets.TSI, 1), size(targets.TSI, 2)), SVector(1, size(targets.TSI, 2))]
    push!(corners, first(corners))
    ij_to_xy(c) = reverse(tform(c))
    # draw extent of TSI within Phantom image
    horiz_mid = [mean((corners[1], corners[2])), mean((corners[3], corners[4]))]
    vert_mid = [mean((corners[1], corners[4])), mean((corners[2], corners[3]))]
    sc = scatter!(ax,ij_to_xy.(corners), color=:red, linewidth=2)
    lines!(ax, ij_to_xy.(corners), color=:red, linewidth=2)
    lines!(ax, ij_to_xy.(horiz_mid), color=:red, linewidth=1)
    lines!(ax, ij_to_xy.(vert_mid), color=:red, linewidth=1)

    axislegend(ax, [sc], ["TSI extent"], position=(:left, :top))

    # set secondary z-axis in inches
    @unpack origin_mm, mm_per_px = calib.phantomscale
    origin_inch = ustrip(u"inch", origin_mm * u"mm")
    inch_per_px = ustrip(u"inch", mm_per_px * u"mm")
    phantom_z_inch = range(origin_inch, length=size(targets.Phantom, 1), step=-inch_per_px)
    phantom_x_inch = range(0, length=size(targets.Phantom, 2), step=inch_per_px)
    zlims_inch = extrema(phantom_z_inch)
    ax_inch = Axis(f[1, 1], yaxisposition=:right, aspect=DataAspect(), 
                yticks=range(round.(Int, zlims_inch)..., step=1), yminorticks=IntervalsBetween(8),
                yminorticksvisible=true, ylabel="z [inch]")
    image!(ax_inch, phantom_x_inch, phantom_z_inch, fill(NaN, size(tgts')))
    hidexdecorations!(ax_inch)

    save(plotsdir("registration", runname(calib) * ".png"), f)
end
# @time Phantom_warp = warp(targets.Phantom, tform, axes(targets.TSI)) |> imadjust
