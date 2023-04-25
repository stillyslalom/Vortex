using DrWatson
@quickactivate "Vortex"
using GLMakie, JLD2, ImageFiltering

# Load PIV runs
runlist = loadmeta(m -> !ismissing(m.TSI_ID) && !ismissing(m.registration_path))

##
function vorticity(xx, yy, u, v)
    Δx, Δy = step(xx), step(yy)
    g₁, g₂ = parent.(Kernel.bickley())
    ImageFiltering.mapwindow(StructArray(; u, v), (3, 3)) do w
        mapreduce(+, eachindex(w)) do i
            w[i].v * g₁[i] / Δx - w[i].u * g₂[i] / Δy
        end
    end
end

runmeta = select_run(runlist, "2023-01-10_CF4_IC2")
PIVcores = map(eachrow(runlist)) do runmeta
    data, path = produce_or_load(runmeta, datadir("PIV", "cores");
                     filename=runname, tag=false) do runmeta
        PIV = load(datadir("PIV", "registered", runname(runmeta)*".jld2"))

        ω = vorticity(PIV["x"], PIV["y"], PIV["u"], PIV["v"])
        ωmax = maximum(abs.(quantile(filter(!isnan, ω), (0.001, 0.999))))

        f = Figure()
        ax = Axis(f[1, 1], aspect=DataAspect())
        p = select_point(ax, marker='●')
        c = [Observable(Point2f(NaN, NaN)), Observable(Point2f(NaN, NaN))]
        heatmap!(ax, PIV["x"], PIV["y"], ω, colormap=:balance, colorrange=(-ωmax, ωmax))
        scatter!(ax, c[1], color=:red)
        scatter!(ax, c[2], color=:blue)

        on(events(f).mousebutton, priority = 100) do event
            if event.button == Mouse.left && event.action == Mouse.press
                if Keyboard.a in events(f).keyboardstate
                    c[1][] = p[]
                    notify(c[1])
                elseif Keyboard.d in events(f).keyboardstate
                    c[2][] = p[]
                    notify(c[2])
                end
            end
            # return Consume(false)
        end

        wait(GLMakie.Screen(f.scene))
        leftcore = (c[1][]...,)
        rightcore = (c[2][]...,)
        xlims = ax.xaxis.attributes.limits[]
        ylims = ax.yaxis.attributes.limits[]
        Dict("leftcore" => leftcore, "rightcore" => rightcore, "xlims" => xlims, "ylims" => ylims)
    end
end

##
PIV = load(datadir("PIV", "registered", runname(runmeta)*".jld2"))

ω = vorticity(PIV["x"], PIV["y"], PIV["u"], PIV["v"])
ωmax = maximum(abs.(quantile(filter(!isnan, ω), (0.001, 0.999))))

f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect())
heatmap!(ax, PIV["x"], PIV["y"], ω, colormap=:balance, colorrange=(-ωmax, ωmax))
scatter!(ax, c[1], color=:red)
scatter!(ax, c[2], color=:blue)

on(events(f).mousebutton, priority = 100) do event
    if event.button == Mouse.left && event.action == Mouse.press
        if Keyboard.a in events(f).keyboardstate
            c[1][] = p[]
            notify(c[1])
        elseif Keyboard.d in events(f).keyboardstate
            c[2][] = p[]
            notify(c[2])
        end
    end
    # return Consume(false)
end

leftcore = (c[1][]...,)
rightcore = (c[2][]...,)
xlims = ax.xaxis.attributes.limits[]
ylims = ax.yaxis.attributes.limits[]
wait(GLMakie.Screen(f.scene))
Dict("leftcore" => leftcore, "rightcore" => rightcore, "xlims" => xlims, "ylims" => ylims)