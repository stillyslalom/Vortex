using DrWatson
@quickactivate "Vortex"
using Vortex, DataFramesMeta, GLMakie, Interpolations
using ImageMorphology, KissSmoothing, NoiseRobustDifferentiation
using ImageFiltering, Statistics, Polynomials
using PyThermo, PyThermo.ShockTube

includet(srcdir("imageutils.jl"))
includet(srcdir("registration.jl"))

savefigs = true

## PLIF background subtraction
PLIFlist = loadmeta() do m
    valid = !ismissing(m.cine_ID)
    valid &= !ismissing(m.timings_path)
end
add_registrations!(PLIFlist)

## Write summary plot
let PLIFmeta = PLIFmeta
# foreach(eachrow(PLIFlist)) do PLIFmeta
    plotpath = plotsdir("mean_spanwise_PLIF", runname(PLIFmeta)*".png")
    # isfile(plotpath) && return nothing
    timings = read_timings(PLIFmeta)
    # time between trigger (if any) and rupture laser pulse
    rupture_lag = (timings[:Rupture1].Shock + timings[:Rupture1].Delay)
    # First Phantom frame occurs at t0_Phantom after rupture
    t0_Phantom = timings[:Death_Star].Delay - rupture_lag

    cine = phantom_bgsub(PLIFmeta)
    I_avg = along(mean, cine, dims=2)
    t_Phantom = range(t0_Phantom + 50, length=size(I_avg, 2), step=50) ./ 1000
    @unpack origin_mm, mm_per_px = PLIFmeta.phantomscale
    z_Phantom = range(stop=origin_mm, step=mm_per_px, length=size(I_avg, 1))

    f = Figure()
    ax1 = Axis(f[1, 1], title=runname(PLIFmeta), xlabel="t [ms]", ylabel="z [mm]")
    hm = heatmap!(ax1, t_Phantom, z_Phantom, rotr90(I_avg))
    savefigs && save(plotpath, f)
    f
end

## Compare shock location between x-t and PLIF
function PLIF_trajectory(PLIFmeta; savefigs=false)
    plotpath = plotsdir("trajectory_PLIF", runname(PLIFmeta)*".png")
    # isfile(plotpath) && return nothing

    # Load x-t data
    xt = load_xt(PLIFmeta)
    shock_x = linear_interpolation(xt.t, xt.x[end] .- xt.x)

    # Load timing information
    timings = read_timings(PLIFmeta)
    shock_times = range(timings[:Death_Star].Delay/1e6, stop=xt.t[end], step=50e-6) # relative to trigger
    shockpos = shock_x.(shock_times)

    # First Phantom frame occurs at t0_Phantom after rupture
    t0_Phantom = timings[:Death_Star].Delay + 50
    t_TSI = (timings[:PIV_trig].Delay + PLIFmeta.PIV_delay) / 1e6

    cine = phantom_bgsub(PLIFmeta)
    I_avg = along(mean, cine, dims=2)
    t_Phantom = range(t0_Phantom, length=size(I_avg, 2), step=50) / 1e6
    @unpack origin_mm, mm_per_px = PLIFmeta.phantomscale
    z_Phantom = range(stop=origin_mm, step=mm_per_px, length=size(I_avg, 1))

    I_median = mapwindow(median, I_avg, (31, 1))
    I_thresh = quantile(vec(I_median[:, eachindex(shockpos)]), 0.1)
    # ringidx = [findfirst(>(max(10/2^12, I_thresh[i])), I_median[:,i-1]) for i in 5:length(shock_times)]
    ringidx = [findfirst(>(max(10/2^12, I_thresh)), I_median[:,i-1]) for i in 5:length(shock_times)]
    ringidx[ringidx .== nothing] .= length(z_Phantom)
    ringpos = reverse(z_Phantom)[ringidx]
    ringpos_smooth, noise = denoise(ringpos; factor = 0.1, rtol=1e-4)

    i_SVI = findmax(ringpos_smooth)[2] + 5 # shock-vortex interaction index
    i_OOW = findlast(ringpos_smooth .> minimum(z_Phantom)) # out-of-window index
    # fit 2nd degree polynomial to the ring position 10 frames before i_SVI
    poly_ring = fit(shock_times[i_SVI-11:i_SVI-1], ringpos_smooth[(i_SVI-11:i_SVI-1) .- 4] .* 1e-3, 1)
    poly_shock = fit(shock_times, shockpos, 1)
    poly_post = fit(shock_times[i_SVI+2:i_OOW], ringpos_smooth[(i_SVI+2:i_OOW) .- 4] .* 1e-3, 1)
    V_pre = Polynomials.derivative(poly_ring)(shock_times[i_SVI-1])
    V_post = Polynomials.derivative(poly_post)(mean(shock_times[i_SVI+2:i_OOW]))
    t_SVI = only(roots(poly_ring - poly_shock))
    dt_TSI = t_TSI - t_SVI

    if savefigs
        f = Figure()
        ax1 = Axis(f[1, 1], title=runname(PLIFmeta), xlabel="t [ms]", ylabel="z [mm]")
        hm = heatmap!(ax1, t_Phantom*1e3, z_Phantom, rotr90(I_avg), colormap=:dense)
        shockplot = lines!(ax1, shock_times[5:end]*1e3, shockpos[5:end]*1e3, color=:red)
        scatter!(ax1, shock_times[5:end]*1e3, ringpos, color=:black, marker=:cross)
        ringplot = lines!(ax1, shock_times[5:end]*1e3, ringpos_smooth, color=:black)
        TSIplot = vlines!(ax1, t_TSI*1e3, color=:green)
        fitplot = lines!(ax1, shock_times[5:end]*1e3, poly_ring.(shock_times[5:end]).*1e3, color=:gray)
        lines!(ax1, shock_times[5:end]*1e3, poly_post.(shock_times[5:end]).*1e3, color=:gray)
        vlines!(ax1, t_SVI*1e3, color=:black)
        xlims!(ax1, (t_Phantom[5], maximum(shock_times)) .* 1e3)
        ylims!(ax1, extrema(z_Phantom))
        axislegend(ax1, [shockplot, ringplot, TSIplot, fitplot], 
                        ["Shock wave", "Vortex ring leading pole", "PIV imaging time", "Least-squares linear fit"],
                         position=(:left, :top))
        save(plotpath, f)
    end
    (; V_pre, V_post, t_SVI, dt_TSI)
end

##
shockruns = filter(r -> !ismissing(r.ptrace_path), PLIFlist)
foreach(eachrow(shockruns)) do PLIFmeta
    PLIF_trajectory(PLIFmeta, savefigs=true)
end
##
# let PLIFmeta = select_run(PLIFlist, "2023-01-18_run8")

transform!(shockruns, AsTable(All()) => ByRow(PLIF_trajectory) => AsTable)
filter!(r -> r.V_post < 0, shockruns)

##
function postshock_u2(PLIFmeta)
    # Load x-t data
    xt = load_xt(PLIFmeta)
    wavespeed = ((xt.x[end] - xt.x[3]) / (xt.t[end] - xt.t[3]))*u"m/s"
    ambient = Species("N2", T = 292, P = 99.629e3)
    shocked, u2 = shockjump(ambient, wavespeed / soundspeed(ambient))
    return ustrip(u"m/s", u2)
end
transform!(shockruns, AsTable(All()) => ByRow(postshock_u2) => :u2)
gb = groupby(shockruns, :MST_gas)
f = Figure()
ax = Axis(f[1, 1], xlabel="Pre-shock velocity [m/s]", ylabel="Post-shock velocity [m/s]")
sc = map(pairs(gb)) do (MST, df)
    scatter!(ax, df.V_pre, df.u2 .+ df.V_post, label=string(MST.MST_gas))
end
axislegend()#, sc, position=(:right, :top))
f
##
PLIFmeta = eachrow(filter(r -> !ismissing(r.ptrace_path), PLIFlist))[11]

timings = read_timings(PLIFmeta)


cine = phantom_bgsub(PLIFmeta)
@unpack origin_mm, mm_per_px = PLIFmeta.phantomscale
x_Phantom = range(0, step=mm_per_px, length=size(cine, 2))
z_Phantom = range(stop=origin_mm, step=mm_per_px, length=size(cine, 1))
f, ax, hm = heatmap(x_Phantom, z_Phantom, rotr90(cine[:,:,48]), interpolate=true, axis=(; aspect = DataAspect()))
ylims!(ax, extrema(z_Phantom))
hlines!(ax, 1e3*shockpos[49], color=:red)
f
##


##
I_total = along(sum, phantom_bgsub(eachrow(PLIFlist)[120]), dims=(1, 2))
I_smooth, noise = denoise(Float64.(I_total); factor = 0.9, rtol=1e-4)
f, ax, sc = scatter(I_total, marker='⋅')
lines!(ax, I_smooth, color=:red)
f

ax2 = Axis(f[2, 1])
lines!(ax2, tvdiff(I_total, 50, 10))
f
##
phantom = phantom_bgsub(eachrow(PLIFlist)[31])
# density(phantom[:,:,10] |> vec, boundary=(0.0, 1))

litreg = along(mean, @views(phantom[:,:,1:150]), dims=3) .> 2*2^-12
closing!(litreg)
# opening!(litreg)
# heatmap(litreg)
# heatmap(along(mean, phantom .* litreg, dims=3))