using DrWatson
@quickactivate "Vortex"
using Vortex
using XLSX, DataFramesMeta
using ImageCore, ImageIO, TiffImages
using ImageTransformations, CoordinateTransformations, StaticArrays
using MATLAB, Statistics, Unitful, GLMakie

includet(srcdir("registration.jl"))
includet(srcdir("imageutils.jl"))

##
function load_targets(tgtdata, tgtdir)
    TSIpath = joinpath(tgtdir, "TSI.png")
    Phantompath = joinpath(tgtdir, "Phantom.png")
    if !isfile(TSIpath)
        TSI = imadjust(rotl90(load(rawdatadir(tgtdata.Date, tgtdata.TSI_tgt_path))))
        save(TSIpath, imadjust(TSI))
    else
        TSI = load(TSIpath)
    end
    if !isfile(Phantompath)
        Phantom = imadjust(Gray.(rotl90(load(rawdatadir(tgtdata.Date, tgtdata.Phantom_tgt_path)))))
        save(Phantompath, imadjust(Phantom))
    else
        Phantom = load(Phantompath)
    end
    return (; TSI, Phantom)
end

function run_cpselect(targets, tgtdir)
    cd(wd) do
        TSIpath = joinpath(tgtdir, "TSI.png")
        Phantompath = joinpath(tgtdir, "Phantom.png")

        mat"Phantom = imread($Phantompath)"
        mat"TSI = imread($TSIpath)"
        mat"[m, f] = cpselect(TSI, Phantom, 'Wait', true)"
        # mat"m = cpcorr(m, f, TSI, Phantom)"
        mat"trans = fitgeotrans(fliplr(m), fliplr(f), 'projective')"
        mat"$tform = trans.T"
        tform
    end
end

## Measure origin, pixel scale, and rotation of Phantom target image
# Select two points along ruler edge of Phantom target image
# and set the distance between them in inches
function selectphantomscale(phantom)
    # Display image with image-standard coordinate axes (origin at top left)
    f, ax, img = image(phantom'; axis=(aspect=DataAspect(), yreversed=true))
    p = select_point(ax, marker='●')
    c = [Observable(p[]), Observable(p[])]
    f[1, 2] = buttongrid = GridLayout(tellheight = false)
    Label(buttongrid[1, 1], "Click to set point selection")
    b  = buttongrid[2:3, 1] = [Button(f, label=@lift(string($(c[1]))), strokecolor=:orange)
                            Button(f, label=@lift(string($(c[2]))), strokecolor=:green)]
    Label(buttongrid[4, 1], "Set point locations on ruler")
    loc_str = buttongrid[5:6, 1] = [Textbox(f, placeholder="Point 1", bordercolor=:orange, validator=Float64),
                                    Textbox(f, placeholder="Point 2", bordercolor=:green, validator=Float64)]
    loc = [Observable(0.0), Observable(0.0)]
    
    for i in 1:2
        on(b[i].clicks) do n
            c[i][] = p[]
            notify(c[i])
        end
        on(loc_str[i].stored_string) do s
            loc[i][] = parse(Float64, s)
        end
    end

    scatter!(ax, p, color=:red, linewidth=2, marker='●')
    scatter!(ax, c[1], color=:orange, markersize=20, marker='+')
    scatter!(ax, c[2], color=:green, markersize=20, marker='+')
    text!(ax, c[1], text=@lift(string('\t', $(loc[1]))), color=:orange, fontsize=20, align=(:right, :center))
    text!(ax, c[2], text=@lift(string('\t', $(loc[2]))), color=:green, fontsize=20, align=(:right, :center))
    wait(GLMakie.Screen(f.scene)) # wait for plot to close

    return (point1 = (coords = c[1][], pos=loc[1][]),
            point2 = (coords = c[2][], pos=loc[2][]))
end

function calcphantomscale(point1, point2)
    px_per_inch = abs(dist2d(point1.coords, point2.coords) / 
                            (point1.pos - point2.pos))*u"inch^-1"
    mm_per_px = ustrip(u"mm", inv(px_per_inch))
    midpoint_px = mean((point1.coords, point2.coords))
    midpoint_mm = ustrip(u"mm", mean((point1.pos, point2.pos)) * u"inch")
    origin_mm = midpoint_mm + midpoint_px[2]*mm_per_px
    return (; origin_mm, mm_per_px)
end

# =============================================================================
# Run calibration for all target images
# =============================================================================
targetlist = DataFrame(XLSX.readtable(datadir("targets.xlsx"), "Sheet1"))

calibration_dicts = map(eachrow(targetlist)) do targetdata
    data, path = produce_or_load(targetdata, datadir("calibrations");
                     filename=runname, tag=false) do targetdata
        calibration_datadir = datadir("calibrations", runname(targetdata))
        # remove & re-make calibration plot directory if it already exists
        isdir(calibration_datadir)  && rm(calibration_datadir; recursive=true)
        mkdir(calibration_datadir)
        targets = load_targets(targetdata, calibration_datadir)

        tf = run_cpselect(targets, calibration_datadir)

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

    f = Figure(resolution=(500, 600))
    ax = Axis(f[1, 1], aspect=DataAspect(), yreversed=true, 
            xlabel="jᵗʰ Phantom pixel", ylabel="iᵗʰ Phantom pixel",
            title = "Phantom and TSI images (warped to Phantom coordinates)\n for $(runname(calib))")
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
