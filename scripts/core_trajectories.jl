using DrWatson
@quickactivate "Vortex"
using Vortex
using GLMakie
using OrderedCollections
using BasicInterpolators
# includet(srcdir("PLIF", "tracking.jl"))
includet(srcdir("imageutils.jl"))

PLIFlist = loadmeta() do m
    valid = !ismissing(m.cine_ID)
    valid &= !ismissing(m.timings_path)
end
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
end
