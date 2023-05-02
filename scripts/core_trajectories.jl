using DrWatson
@quickactivate "Vortex"
using Vortex
using GLMakie
using OrderedCollections
using BasicInterpolators
using Dierckx
using LsqFit
using PyThermo
# includet(srcdir("PLIF", "tracking.jl"))
includet(srcdir("imageutils.jl"))

PLIFlist = loadmeta() do m
    valid = !ismissing(m.cine_ID)
    valid &= !ismissing(m.timings_path)
    valid &= !ismissing(m.registration_path)
end
Vortex.add_registrations!(PLIFlist)
##
meta = select_run(PLIFlist, "2023-01-23_run3")

PLIFcores = map(eachrow(PLIFlist)) do runmeta
    data, path = produce_or_load(runmeta, datadir("PLIF", "cores");
    filename=runname, tag=false) do runmeta
        f = Figure(resolution=(1100,1300))
        ax = Axis(f[1, 1], aspect=DataAspect(), title=runname(runmeta))
        PLIF = eachslice(phantom_bgsub(runmeta), dims=3)
        sg = SliderGrid(f[2, 1], (label = "Index", range=eachindex(PLIF), startvalue=1))
        frame = @lift(rotr90(PLIF[$(sg.sliders[1].value)]))
        vbounds = @lift(quantile($frame, (0.001, 0.9999)))
        heatmap!(ax, frame, colorrange=vbounds, colormap=:grays)

        # Core position storage
        core = OrderedDict(
            :left => OrderedDict(
                :preshock =>    OrderedDict{Int,Point2f}(), 
                :postshock =>   OrderedDict{Int,Point2f}(), 
                :postreshock => OrderedDict{Int,Point2f}()),
            :right => OrderedDict(
                :preshock =>    OrderedDict{Int,Point2f}(), 
                :postshock =>   OrderedDict{Int,Point2f}(), 
                :postreshock => OrderedDict{Int,Point2f}()))

        lcore_menu = Menu(f, options=zip(["Pre-shock", "Post-shock", "Post-reshock"], core[:left]), width=120)
        rcore_menu = Menu(f, options=zip(["Pre-shock", "Post-shock", "Post-reshock"], core[:right]), width=120)
        lcore_delete = Button(f, label="Delete", halign=:left)
        rcore_delete = Button(f, label="Delete", halign=:left)
        zoom_menu = Toggle(f)
        f[1, 2] = vgrid!(Label(f, "Left core", width=nothing), lcore_menu, lcore_delete,
                        Label(f, "Right core", width=nothing), rcore_menu, rcore_delete;
                        tellheight=false, width=180)
        f[2, 2] = hgrid!(Label(f, "Zoom", width=nothing), zoom_menu; tellheight=false, width=110)

        lcores = @lift(last($(lcore_menu.selection)))
        rcores = @lift(last($(rcore_menu.selection)))

        lcore_sc = scatter!(ax, @lift(collect(values($lcores))), color=:green, markersize=10)
        rcore_sc = scatter!(ax, @lift(collect(values($rcores))), color=:red, markersize=10)

        # Interpolated track plotting
        lines!(ax, @lift(collect(values($lcores))), color=:green, markersize=10)
        lines!(ax, @lift(collect(values($rcores))), color=:red, markersize=10)

        scatter!(ax, @lift(get($lcores, $(sg.sliders[1].value), Point2f(NaN, NaN))), 
            color=:green, markersize=30, marker='+')
        scatter!(ax, @lift(get($rcores, $(sg.sliders[1].value), Point2f(NaN, NaN))),
            color=:red, markersize=30, marker='+')

        p = select_point(ax, marker='+', markersize=36)
        # Keyboard controls
        on(events(f).keyboardbutton) do event
            if event.action == Keyboard.press || event.action == Keyboard.repeat
                if (event.key == Keyboard.right) || (event.key == Keyboard.d) 
                    set_close_to!(sg.sliders[1], sg.sliders[1].value[] + 1)
                elseif (event.key == Keyboard.left) || (event.key == Keyboard.a)
                    set_close_to!(sg.sliders[1], sg.sliders[1].value[] - 1)
                elseif (event.key == Keyboard.l) || (event.key == Keyboard.q)
                    ldict = lcores.val
                    ldict[sg.sliders[1].value[]] = p[]
                    sort!(ldict)
                    notify(lcore_menu.selection)
                elseif event.key == Keyboard.e
                    rdict = rcores.val
                    rdict[sg.sliders[1].value[]] = p[]
                    sort!(rdict)
                    notify(rcore_menu.selection)
                end
            end
        end

        # delete core
        on(lcore_delete.clicks) do c
            delete!(lcores[], sg.sliders[1].value[])
            notify(lcore_menu.selection)
        end

        on(rcore_delete.clicks) do c
            delete!(rcores[], sg.sliders[1].value[])
            notify(rcore_menu.selection)
        end

        # Enable/disable drag-to-zoom
        on(zoom_menu.active) do active
            if active[]
                activate_interaction!(ax, :rectanglezoom)
            else
                deactivate_interaction!(ax, :rectanglezoom)
            end
        end
        notify(zoom_menu.active)

        wait(GLMakie.Screen(f.scene))
        return core
    end
    data["runname"] = runname(runmeta)
    data
end

##
function recursive_dict_count(d)
    # count number of elements in a nested dictionary
    if isa(d, AbstractDict)
        return sum(map(recursive_dict_count, values(d)))
    else
        return 1
    end
end

mapreduce(recursive_dict_count, +, PLIFcores) # 4805 core positions

## PLIFcores entry cleaning
# Subsequent cores with identical values and different keys are invalid and should be removed
function clean_cores!(cores)
    for side in keys(cores)
        for state in keys(cores[side])
            prevval = Point2f(NaN, NaN)
            for (i, (k, v)) in enumerate(cores[side][state])
                if v == prevval
                    println("deleting $side $state $k")
                    delete!(cores[side][state], k)
                end
                prevval = v
            end
        end
    end
    cores
end

# test data
testcores = Dict(
    "left" => OrderedDict(
        :preshock => OrderedDict(
            1 => Point2f(1, 1),
            2 => Point2f(1, 1), # duplicate
            3 => Point2f(3, 3),
        )
    )
)
recursive_dict_count(testcores)
clean_cores!(testcores) |> recursive_dict_count

##
PLIFlist.runname = runname.(eachrow(PLIFlist))
leftjoin!(PLIFlist, DataFrame(PLIFcores), on=:runname)

##
function MST_state(runmeta)
    p_ace = Species("acetone", T = 22u"°C").Psat*u"Pa"
    p_total = (runmeta.MST_psig + 14.5)*u"psi" |> u"Pa"
    χ_ace = p_ace/p_total
    MSTgas = Mixture([runmeta.MST_gas => 1-χ_ace, "acetone" => χ_ace])
    Ma = fzero(1.31) do M
        pressure(PyThermo.ShockTube.driverpressure(MSTgas, MSTgas, M)) - p_total
    end
    sc = PyThermo.ShockTube.shockcalc(MSTgas, MSTgas, Ma)
end

# runmeta = eachrow(PLIFlist)[5]
# begin; runmeta = select_run(PLIFlist, "2022-11-15_run1")
savefigs = true
foreach(eachrow(PLIFlist)) do runmeta
    cine = Vortex.phantom_bgsub(runmeta)
    i_PIV, t_err = Vortex.phantom_timing(runmeta, Vortex.read_timings(runmeta), size(cine, 3))
    t_err > 0 && @warn "PIV time is not simultaneous with cine frame"
    cine_PIV = cine[:, :, i_PIV]
    phantomscale = runmeta.phantomscale
    xx = range(0, step=phantomscale.mm_per_px, length=size(cine_PIV, 2)) ./ 1000
    y0 = 0.01425 # distance from end wall of shock tube to tip of vortex generator
    yy = range(stop = phantomscale.origin_mm, step=phantomscale.mm_per_px, length=size(cine_PIV, 1)) ./ 1000 .- y0

    sc = MST_state(runmeta)
    Ws = soundspeed(sc.driven)*sc.Ms # shock speed in driven gas
    t_shockprop = 10.61u"inch" / Ws

    timings = read_timings(runmeta)
    rupture_lag = (timings[:Rupture1].Shock + timings[:Rupture1].Delay)*u"μs" + t_shockprop
    # First Phantom frame occurs at t0_Phantom after rupture
    t0_Phantom = timings[:Death_Star].Delay - round(Int, ustrip(u"μs", rupture_lag))

    t_Phantom = range(t0_Phantom + 50, length=size(cine, 3), step=50) ./ 1e6

    itp_xgrid = Spline1D(eachindex(xx), xx, k=1, bc="extrapolate")
    itp_ygrid = Spline1D(eachindex(yy), yy, k=1, bc="extrapolate")

    PLIFcore = runmeta
    core_ix(PLIFcore, lr, state) = first.(values(PLIFcore[lr][state]))
    core_iy(PLIFcore, lr, state) = last.(values(PLIFcore[lr][state]))
    core_it(PLIFcore, lr, state) = collect(keys(PLIFcore[lr][state]))

    f = Figure(resolution = (400, 300), fontsize=12)
    ax = Axis(f[1, 1], title=runname(runmeta) * " -  $(runmeta.MST_gas) at $(runmeta.MST_psig) psig",
        xlabel="time (s)", ylabel="y (m)")

    yfits = Vector{Float64}[]

    isempty(PLIFcore[:left][:preshock]) && isempty(PLIFcore[:right][:preshock]) && return nothing
    for lr in (:left, :right)
        t = t_Phantom[core_it(PLIFcore, lr, :preshock)]
        isempty(t) && continue
        x = itp_xgrid(core_ix(PLIFcore, lr, :preshock))
        y = itp_ygrid(core_iy(PLIFcore, lr, :preshock))

        @. y_c(i, p) = p[1]*(1 - exp(-p[2]*(i - p[3])))
        yfit = curve_fit(y_c, [0; t], [0; y], [0.5, 200, 1e-3])
        # yfit = curve_fit(y_c, t, y, [0.5, 200, 1e-3])
        push!(yfits, yfit.param)

        scatter!(ax, t, y, label="Pre-shock " * string(lr))
        tfine = LinRange(0, maximum(t), 100)
        lines!(ax, tfine, y_c(tfine, yfit.param))
    end

    psfits = Vector{Float64}[]
    if !isempty(PLIFcore[:left][:postshock]) || !isempty(PLIFcore[:right][:postshock])
        for lr in (:left, :right)
            x = itp_xgrid(core_ix(PLIFcore, lr, :postshock))
            y = itp_ygrid(core_iy(PLIFcore, lr, :postshock))
            t = t_Phantom[core_it(PLIFcore, lr, :postshock)]

            @. y_c(i, p) = p[1] + p[2]*i
            yfit = curve_fit(y_c, t, y, Float64[1, -200])
            # yfit = curve_fit(y_c, t, y, [0.5, 200, 1e-3])
            push!(psfits, yfit.param)

            scatter!(ax, t, y, label="Post-shock " * string(lr))
            tfine = LinRange(extrema(t)..., 100)
            lines!(ax, tfine, y_c(tfine, yfit.param))
        end
    end
    u₀ = isempty(yfits) ? NaN : mean(yfit -> yfit[1]*yfit[2], yfits)
    u = isempty(yfits) ? NaN : mean(yfit -> yfit[1]*yfit[2]*exp(-yfit[2]), yfits)
    uₚₛ = isempty(psfits) ? NaN : mean(yfit -> yfit[2], psfits)
    xlims!(ax, 0, nothing)
    ylims!(ax, 0, nothing)
    text!(ax, 0.002, 0.1, fontsize=12, text = """
    u₀ = $(round(u₀, digits=1)) m s^-1
    uₚ = $(round(u"m/s", sc.u2, digits=1))
    u₀/uₚ = $(round(u₀*u"m/s" / sc.u2, digits=2))
    """ * (isempty(psfits) ? "" : "uₚₛ = $(round(uₚₛ, digits=1)) m s^-1"))
    try
        axislegend(ax, position=:lt, rowgap=-5)
    catch
        nothing
    end
    save(plotsdir("PLIF_core_trajectories", runname(runmeta)*".png"), f, px_per_unit=4)
    f
end

#
# f = Figure(resolution=(1100,1300))
# ax = Axis(f[1, 1], aspect=DataAspect(), title=runname(runmeta))
# PLIF = eachslice(cine, dims=3)
# sg = SliderGrid(f[2, 1], (label = "Index", range=eachindex(PLIF), startvalue=1))
# frame = @lift(rotr90(PLIF[$(sg.sliders[1].value)]))
# vbounds = @lift(quantile($frame, (0.001, 0.9999)))
# heatmap!(ax, xx, yy, frame, colorrange=vbounds, colormap=:grays)
# f
#