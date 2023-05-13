using DrWatson
@quickactivate "Vortex"
using Vortex
using GLMakie
using OrderedCollections
using BasicInterpolators
using Dierckx
using LsqFit
using PyThermo
using LaTeXStrings
using DataFramesMeta
using Interpolations
using ImageFiltering

# includet(srcdir("PLIF", "tracking.jl"))
includet(srcdir("imageutils.jl"))

core_ix(PLIFcore, lr, state) = first.(values(PLIFcore[lr][state]))
core_iy(PLIFcore, lr, state) = last.(values(PLIFcore[lr][state]))
core_it(PLIFcore, lr, state) = collect(keys(PLIFcore[lr][state]))
##

PLIFlist = loadmeta() do m
    valid = !ismissing(m.cine_ID)
    valid &= !ismissing(m.timings_path)
    valid &= !ismissing(m.registration_path)
end
Vortex.add_registrations!(PLIFlist)

masters = map(eachrow(PLIFlist)) do runmeta
    data = JLD2.load(datadir("master", runname(runmeta)*".jld2"))
    data["runname"] = runname(runmeta)
    data
end |> DataFrame

PLIFlist.runname = runname.(eachrow(PLIFlist))
leftjoin!(PLIFlist, masters, on = :runname)
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

# runmeta = eachrow(PLIFlist)[5]
begin; runmeta = select_run(PLIFlist, "2022-11-15_run1")
savefigs = false
# foreach(eachrow(PLIFlist)) do runmeta
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
    f = Figure(resolution = (400, 300), fontsize=12)
    ax = Axis(f[1, 1], title=runname(runmeta) * " -  $(runmeta.MST_gas) at $(runmeta.MST_psig) psig",
        xlabel=L"$t$ [s]", ylabel=L"$z$ [m]")

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
    u = isempty(yfits) ? NaN : mean(yfit -> yfit[1]*yfit[2]*exp(-yfit[2]*runmeta.t_SVI), yfits)
    uₚₛ = isempty(psfits) ? NaN : mean(yfit -> yfit[2], psfits)
    xlims!(ax, 0, nothing)
    ylims!(ax, 0, nothing)
    texts = [
        L"""$u_0 = %$(round(u₀, digits=1))$ m/s""",
        L"""$u_p = %$(ustrip(round(u"m/s", sc.u2, digits=1)))$ m/s""",
        L"""$u_0/u_p = %$(round(u₀*u"m/s" / sc.u2, digits=2))$"""]
    isempty(psfits) || push!(texts, L"$u^{-} = %$(round(u, digits=1))$ m/s")
    isempty(psfits) || push!(texts, L"$u^{+} = %$(round(uₚₛ, digits=1))$ m/s")
    text!(ax, fill(0.5, length(texts)), reverse(range(0.02, step=0.07, length=length(texts))),
        text = texts, space=:relative, fontsize=12)
    try
        axislegend(ax, position=:lt, rowgap=-5)
    catch
        nothing
    end
    savefigs && save(plotsdir("PLIF_core_trajectories", "z", runname(runmeta)*".png"), f, px_per_unit=4)
    f
end

## Calculate and plot mean ± std PLIF diameter versus time for each group of runs

# Find minimum & maximum time for each core group
function core_timebounds(cores, t0=0.0)
    tmin = mapreduce(c -> c.t - t0, min, values(cores), init=Inf)
    tmax = mapreduce(c -> c.t - t0, max, values(cores), init=-Inf)
    tmin, tmax
end

function bothcore_timebounds(cores, state, t0=0.0)
    ltmin, ltmax = core_timebounds(cores[:left][state], t0)
    rtmin, rtmax = core_timebounds(cores[:right][state], t0)
    min(ltmin, rtmin), max(ltmax, rtmax)
end

transform!(PLIFlist, :cores => ByRow(c -> bothcore_timebounds(c, :preshock)) => :preshock_timebounds)
transform!(PLIFlist, [:cores, :t_SVI] => ByRow((c, t0) -> bothcore_timebounds(c, :postshock, t0)) => :postshock_timebounds)

# Calculate linear interpolants
function xz_linear_itps(cores, xz, t0=0.0)
    T_itp = typeof(linear_interpolation([1, 2.2, 3.2], rand(3), extrapolation_bc=NaN))
    itps = OrderedDict{Tuple{Symbol,Symbol}, T_itp}()
    for lr in (:left, :right)
        for state in (:preshock, :postshock)
            cdict = cores[lr][state]
            t = getfield.(values(cdict), :t)
            x = getfield.(values(cdict), xz)
            if length(t) < 2
                itps[(lr, state)] = LinearInterpolation([0, 0 + eps()], [0, 0 + eps()], extrapolation_bc=NaN)
                continue
            end

            state == :postshock && (t .-= t0)
            itps[(lr, state)] = linear_interpolation(t, x, extrapolation_bc=NaN)
        end
    end
    itps
end

transform!(PLIFlist, [:cores, :t_SVI] => ByRow((c, t0) -> xz_linear_itps(c, :x, t0)) => :x_itps)
transform!(PLIFlist, [:cores, :t_SVI] => ByRow((c, t0) -> xz_linear_itps(c, :z, t0)) => :z_itps)

## find 3rd largest & smallest min/max time bounds for each group
# and set time range accordingly
function kth_timebounds(timebounds)
    tmins = first.(timebounds)
    tmaxs = last.(timebounds)
    sort!(tmins)
    sort!(tmaxs)
    tmins[3], tmaxs[end-2]
end

PLIFlist.shockrun = @. !ismissing(PLIFlist.t_SVI)

gbIC = groupby(PLIFlist, [:MST_gas, :MST_psig])
transform!(gbIC, :preshock_timebounds => kth_timebounds => :preshock_group_timebounds)

PLIFlist.t_IC = [range(t..., step=50e-6) for t in PLIFlist.preshock_group_timebounds]

## Sample diameter at each time point
function ring_diameter(xi, zi, t, state=:preshock)
    xl = xi[(:left, state)](t)
    xr = xi[(:right, state)](t)
    zl = zi[(:left, state)](t)
    zr = zi[(:right, state)](t)
    hypot(xr - xl, zr - zl)
end

## Sample centroid at each time point
function ring_centroid(xi, zi, t, state=:preshock)
    xl = xi[(:left, state)](t)
    xr = xi[(:right, state)](t)
    zl = zi[(:left, state)](t)
    zr = zi[(:right, state)](t)
    (; x = (xl + xr)/2, z = (zl + zr)/2)
end

transform!(PLIFlist, 
    [:t_IC, :x_itps, :z_itps] => ByRow((t, xi, zi) -> ring_diameter.(Ref(xi), Ref(zi), t)) => :D_IC)
transform!(PLIFlist, 
    [:t_IC, :x_itps, :z_itps] => ByRow((t, xi, zi) -> ring_centroid.(Ref(xi), Ref(zi), t)) => :ring_centroid)

## Calculate mean & std diameter at each time point for each group
function mean_std_diameter(D, i)
    finiteD = filter(isfinite, getindex.(D, i))
    μ = mean(finiteD)
    σ = std(finiteD, mean=μ)
    μ, σ
end

function group_mean_std_diameter(D, t)
    map(enumerate(t[1])) do (i, t)
        Ds = getindex.(D, i)
        finiteD = filter(isfinite, Ds)
        μ = mean(finiteD)
        σ = std(finiteD, mean=μ)
        (; t, μ, σ)
    end
end

## Plot mean ± std diameter versus time for each group
IC_D_grps = groupby(combine(gbIC, [:D_IC, :t_IC] => group_mean_std_diameter => :D_IC_mean_std),
    [:MST_gas, :MST_psig])

f = Figure(resolution=(700,600))
gas_locs = Dict("N2" => (1, 1), "Ar" => (1, 2), "CF4" => (2, 1), "SF6" => (2, 2))
gas_clean = Dict("N2" => L"\mathrm{N_2}", "Ar" => L"\mathrm{Ar}", "CF4" => L"\mathrm{CF_4}", "SF6" => L"\mathrm{SF_6}")
axs = [Axis(f[i, j], xlabel=L"$t$ [ms]$", ylabel=L"$D$ [m]") for i in 1:2, j in 1:2]
colorset = Dict(unique(PLIFlist.MST_psig) .=> Makie.wong_colors()[1:7])
lineplots = Dict()
for k in sort(keys(IC_D_grps), by=(k -> k.MST_psig))
    t = getfield.(IC_D_grps[k].D_IC_mean_std, :t)
    t .*= 1e3
    μ = getfield.(IC_D_grps[k].D_IC_mean_std, :μ)
    t[.!isfinite.(μ)] .= NaN
    # μ = imfilter(μ, KernelFactors.gaussian(3))
    σ = clamp01nan!(getfield.(IC_D_grps[k].D_IC_mean_std, :σ))
    # clamp!(σ, 1e-3, )
    # σ = imfilter(σ, KernelFactors.gaussian(3))
    length(t) < 20 && continue
    axloc = gas_locs[k.MST_gas]
    band!(axs[axloc...], t, μ .- σ, μ .+ σ, color=RGBA(RGB(colorset[k.MST_psig]), 0.3))
    p = lines!(axs[axloc...], t, μ, 
        label=L"%$(round((k.MST_psig + 14.5)/14.5, digits=2))",
        color=colorset[k.MST_psig], linewidth=2)
end
for (gas, loc) in gas_locs
    current_axis!(f, axs[loc...])
    axislegend(latexstring(gas_clean[gas], " ", L"p_4/p_1"), rowgap=0)
end
xlims!(axs[1, 1], nothing, 9.5)
xlims!(axs[2, 1], nothing, 24)
xlims!(axs[2, 2], nothing, 24)

save(plotsdir("PLIF_core_trajectories", "core_diameter.svg"), f)

f
## Repeat same plot in normalized dimensions (D/D0), t*u_p/D0
u_p = Dict(k => MST_state(k.MST_gas, k.MST_psig).u2 for k in keys(IC_D_grps))

f = Figure(resolution=(700,600))
axs = [Axis(f[i, j], xlabel=L"$t\, u_p / D_p$", ylabel=L"$D_\mathrm{ring}/D_p$") for i in 1:2, j in 1:2]
lineplots = Dict()
for k in sort(keys(IC_D_grps), by=(k -> k.MST_psig))
    t = getfield.(IC_D_grps[k].D_IC_mean_std, :t)
    t .*= u"s" * u_p[k]/(0.875u"inch" |> u"m")
    μ = getfield.(IC_D_grps[k].D_IC_mean_std, :μ)
    μ .*= u"m" / (0.875u"inch" |> u"m")
    t[.!isfinite.(μ)] .= NaN
    σ = clamp01nan!(getfield.(IC_D_grps[k].D_IC_mean_std, :σ))
    σ .*= u"m" / (0.875u"inch" |> u"m")
    length(t) < 20 && continue
    axloc = gas_locs[k.MST_gas]
    band!(axs[axloc...], t, μ .- σ, μ .+ σ, color=RGBA(RGB(colorset[k.MST_psig]), 0.3))
    p = lines!(axs[axloc...], t, μ, 
        label=L"%$(round((k.MST_psig + 14.5)/14.5, digits=2))",
        color=colorset[k.MST_psig], linewidth=2)
end
for (gas, loc) in gas_locs
    current_axis!(f, axs[loc...])
    axislegend(latexstring(gas_clean[gas], " ", L"p_4/p_1"), rowgap=0, position=:lt)
    xlims!(axs[loc...], -3, nothing)
end
save(plotsdir("PLIF_core_trajectories", "normalized_core_diameter.svg"), f)

f

## Calculate post-shock diameter normalized by pre-shock diameter
gbPS = groupby(filter(r -> r.shockrun, PLIFlist), [:MST_gas, :MST_psig])
PSlist = filter(r -> r.nrow >= 3, transform!(gbPS, nrow))
gbPS = groupby(PSlist, [:MST_gas, :MST_psig])
transform!(gbPS, :postshock_timebounds => kth_timebounds => :postshock_group_timebounds)
PSlist.t_PS = [range(t..., step=5e-6) for t in PSlist.postshock_group_timebounds]

transform!(PSlist, 
    [:t_PS, :x_itps, :z_itps] => ByRow((t, xi, zi) -> ring_diameter.(Ref(xi), Ref(zi), t, :postshock)) => :D_PS)

PSlist.D⁻ = [D[findlast(!isnan, D)] for D in PSlist.D_IC]
PSlist.D̄ = [p.D_PS ./ p.D⁻ for p in eachrow(PSlist)]

PS_D_grps = groupby(combine(gbPS, [:D̄, :t_PS] => group_mean_std_diameter => :D_PS_mean_std),
    [:MST_gas, :MST_psig])

## Plot post-shock diameter versus post-shock time
f = Figure(resolution=(700,500))
ax = Axis(f[1, 1], xlabel=L"$t - t_{SVI}$ [ms]", ylabel=L"$D/D^-$")
lineplots = Dict()
i = 1
gasorder = Dict(["N2", "Ar", "CF4", "SF6"] .=> 1:4)
for k in sort(keys(PS_D_grps), by=(k -> gasorder[k.MST_gas]))
    t = getfield.(PS_D_grps[k].D_PS_mean_std, :t)
    t .*= 1e3
    μ = getfield.(PS_D_grps[k].D_PS_mean_std, :μ)
    ig = isfinite.(μ)
    @show k, any(isnan, μ)
    σ = clamp01nan!(getfield.(PS_D_grps[k].D_PS_mean_std, :σ))
    t = t[ig]
    μ = μ[ig]
    σ = σ[ig]
    length(t) < 20 && continue
    k.MST_psig > 30 && continue
    axloc = gas_locs[k.MST_gas]
    band!(ax, t, μ .- σ, μ .+ σ, color=RGBA(RGB(Makie.wong_colors()[i]), 0.3))
    p = lines!(ax, t, μ, 
        label=gas_clean[k.MST_gas],
        color=Makie.wong_colors()[i], linewidth=2)
    i += 1
end
axislegend(ax, rowgap=0, position=:lt)

save(plotsdir("PLIF_core_trajectories", "postshock_diameter.svg"), f)
f
