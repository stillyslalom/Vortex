using ImageIO, TiffImages, MATLAB, ImageCore, Statistics, GLMakie, Unitful

include(srcdir("imageutils.jl"))

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

"Distance between two points"
dist2d(x, y) = sqrt(sum((x .- y).^2))
dist2d(x) = dist2d(x...)

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

