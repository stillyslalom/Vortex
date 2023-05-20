using DrWatson
@quickactivate "Vortex"
using Vortex
using GLMakie, CairoMakie
using JLD2
using Measurements
using LsqFit
using LaTeXStrings
using Printf
using Unitful
using PyThermo

update_theme!(Theme(fonts = (; regular = "New Roman", bold = "New Roman Bold")))
##
corefield(cores, lr, state, field) = getfield.(values(cores[lr][state]), field)
##

runlist = loadmeta() do m
    !ismissing(m.TSI_ID) && !ismissing(m.registration_path) && !ismissing(m.cine_ID) && !ismissing(m.timings_path)
end

masters = map(eachrow(runlist)) do runmeta
    data = JLD2.load(datadir("master", runname(runmeta)*".jld2"))
    data["runname"] = runname(runmeta)
    data
end |> DataFrame

runlist.runname = runname.(eachrow(runlist))
leftjoin!(runlist, masters, on = :runname)

##
runlist.At = map(eachrow(runlist)) do s
    ρvortex = PyThermo.density(MST_state(s).driven)
    ρambient = PyThermo.density(Species("N2", P=14.3u"psi", T=19u"°C"))
    At = (ρvortex .- ρambient) ./ (ρvortex .+ ρambient)
end

##
# runmeta = select_run(runlist, "2023-02-10_Ar_IC9")
# PIV = load(datadir("PIV", "registered", runname(runmeta)*".jld2"))
# ω = vorticity(PIV["u"], PIV["v"], step(runmeta.grid.x), step(runmeta.grid.z))
# heatmap(ω, axis=(; aspect=DataAspect()), colormap=:balance, colorrange=(-1e5, 1e5))

##
# ω = vorticity(PIV["u"] .± PIV["σu"], PIV["v"] .± PIV["σv"], step(runmeta.grid.x), step(runmeta.grid.z))
∫s = map(eachrow(runlist)) do runmeta
    corepath = datadir("PIV", "cores", runname(runmeta)*".jld2")
    !isfile(corepath) && return missing
    core = load(corepath)
    isnan(core["leftcore"][1]) || isnan(core["rightcore"][1]) && return missing
    centerpoint = (core["leftcore"] .+ core["rightcore"]) ./ 2

    PIV = load(datadir("PIV", "registered", runname(runmeta)*".jld2"))
    ω = Vortex.vorticity(PIV["u"], PIV["v"], step(runmeta.grid.x), step(runmeta.grid.z))

    i_x = findfirst(x -> x > core["xlims"][1], PIV["x"]):findlast(x -> x < core["xlims"][2], PIV["x"])
    i_y = findfirst(y -> y > core["ylims"][1], PIV["y"]):findlast(y -> y < core["ylims"][2], PIV["y"])

    xx = PIV["x"][i_x] .* 1e-3 .* u"m"
    yy = PIV["y"][i_y] .* 1e-3 .* u"m"
    ωᵥ = ω[i_x, i_y] .* u"s^-1"
    ωᵥ[isnan.(ωᵥ)] .= 0u"s^-1"
    clamp!(ωᵥ, quantile(vec(ωᵥ), (0.0001, 0.9999))...)

    function vortex_integrals(ω, r, z)
        dr = step(r)
        dz = step(z)
        dΓ = ω .* dr * dz
        Γ = sum(dΓ)
        R = sum(dΓ .* r) / Γ
        Z = sum(dΓ .* z') / Γ
        R₂² = sum(dΓ .* r.^2) / Γ
        R_diff = R₂² - R^2
        a = R_diff < 0u"m^2" ? R/2 : sqrt(R_diff)
        (; Γ, R, Z, R₂², a)
    end

    _, i_xmid = findmin(x -> abs(x - centerpoint[1]*u"mm"), xx)
    xmid = xx[i_xmid]

    ∫ᵣ = vortex_integrals(ωᵥ[i_xmid:end, :], xx[i_xmid:end] .- xmid, yy)
    ∫ₗ = vortex_integrals(ωᵥ[1:i_xmid, :], xx[1:i_xmid] .- xmid, yy)
    Γ = (∫ₗ.Γ - ∫ᵣ.Γ)/2
    (; runname = runname(runmeta), ∫ᵣ, ∫ₗ, Γ)
end |> DataFrame

leftjoin!(runlist, ∫s, on = :runname)

##
runlist.PIV_quality = map(eachrow(runlist)) do runmeta
    maskpath = datadir("PIV", "masks", runname(runmeta)*".jld2")
    !isfile(maskpath) && return missing
    load(maskpath, "quality")
end

##
good_PIV = filter(r -> r.PIV_quality > 0, runlist)
good_PIV.shockrun = .!ismissing.(good_PIV.ptrace_path)
good_PIV.shockrun .&= ifelse.(good_PIV.shockrun, good_PIV.t_SVI .< good_PIV.t_TSI, false)

gdf = groupby(good_PIV, [:MST_gas, :MST_psig, :shockrun])
shockICs = ["N2" => 20.5, "Ar" => 20.5, "CF4" => 25.5, "SF6" => 25.5]

@. Γ_model(t, p) = p[1] + t * p[2]
Γ_fits = Dict{String,Vector{Float64}}()

f = Figure(resolution=100 .* (7, 5))
for j in 1:2, i in 1:2
    local ax = Axis(f[i, j], xlabel="time [ms]", ylabel="ΔΓ [m²/s]", 
        xminorticksvisible=true, yminorticksvisible=true)
    gas, psig = shockICs[j + 2*(i-1)]
    ICs = gdf[(MST_gas = gas, MST_psig=psig, shockrun=false)]
    shocks = gdf[(MST_gas = gas, MST_psig=psig, shockrun=true)]
    t_init = 1e3*ICs.t_TSI
    Γ_init = ustrip.(u"m^2/s", ICs.Γ)
    cf = curve_fit(Γ_model, t_init, Γ_init, ones(2))
    Γ_fits[gas] = cf.param
    lines!(ax, 1e3*ICs.t_TSI, Γ_model.(t_init, Ref(cf.param)), color=:black, label="IC")
    plot!(ax, 1e3*ICs.t_TSI, ustrip.(u"m^2/s", ICs.Γ), color=:black, label="IC", marker=:diamond)
    plot!(ax, 1e3*shocks.t_TSI, ustrip.(u"m^2/s", shocks.Γ), color=:red, marker=:dtriangle, label="shock")
    axislegend(gas, position = i == 1 ? :lb : :lt, merge=true)
end
N2_ax, CF4_ax, Ar_ax, SF6_ax = f.content[1:2:end]
linkyaxes!(N2_ax, Ar_ax)
linkyaxes!(CF4_ax, SF6_ax)
linkxaxes!(N2_ax, CF4_ax)
linkxaxes!(Ar_ax, SF6_ax)
hidexdecorations!.((N2_ax, Ar_ax), ticks=false, grid=false, minorticks=false)
hideydecorations!.((Ar_ax, SF6_ax), ticks=false, grid=false, minorticks=false)

save(plotsdir("circulation_jump", "Gamma_IC_shock.svg"), f,)
f
## Plot circulation jump (relative to pre-shock trend) versus time relative to t_SVI
f = Figure(resolution=100 .* (7, 5))
ax1 = Axis(f[1, 1], xlabel="post-shock time [ms]", ylabel="ΔΓ [m²/s]", 
        xminorticksvisible=true, yminorticksvisible=true)
ΔΓtrend = Dict{String,Vector{Float64}}()
for i in 1:2, j in 1:2
    idx = j + 2*(i-1)
    gas, psig = shockICs[idx]
    # local ax = Axis(f[i, j], xlabel="time [ms]", ylabel="ΔΓ [m²/s]", 
    #     xminorticksvisible=true, yminorticksvisible=true)
    shocks = gdf[(MST_gas = gas, MST_psig=psig, shockrun=true)]
    cf = Γ_fits[gas]
    ΔΓ = shocks.Γ .- Γ_model.(1e3*shocks.t_TSI, Ref(cf))*u"m^2/s"
    ΔΓtrend[gas] = ustrip.(u"m^2/s", ΔΓ)
    plot!(ax1, 1e3*(shocks.t_TSI .- shocks.t_SVI), ustrip.(u"m^2/s", ΔΓ), 
        marker=:dtriangle, label=gas, color=Makie.wong_colors()[idx])
    hlines!(ax1, ustrip(u"m^2/s", mean(ΔΓ)), color=Makie.wong_colors()[idx], linestyle=:dot, label=gas)
end
axislegend(position = :rt, merge=true)
save(plotsdir("circulation_jump", "Delta_Gamma_versus_IC_trend.svg"), f,)
# Plot ΔΓ distribution versus Atwood number
# f = Figure(resolution=100 .* (7, 5))
ax2 = Axis(f[1, 2], xlabel="Atwood number", 
        xminorticksvisible=true, yminorticksvisible=true)
for i in 1:2, j in 1:2
    idx = j + 2*(i-1)
    gas, psig = shockICs[idx]
    ΔΓ = ΔΓtrend[gas]
    shocks = gdf[(MST_gas = gas, MST_psig=psig, shockrun=true)]
    ρvortex = [PyThermo.density(MST_state(s).driven) for s in eachrow(shocks)]
    ρambient = PyThermo.density(Species("N2", P=14.3u"psi", T=19u"°C"))
    At = (ρvortex .- ρambient) ./ (ρvortex .+ ρambient)
    boxplot!(ax2, fill(mean(At), length(ΔΓ)), ΔΓ, 
        label=gas, color=Makie.wong_colors()[idx], width=0.07)
end
linkyaxes!(ax1, ax2)
save(plotsdir("circulation_jump", "Delta_Gamma_versus_At.svg"), f,)
f

## Find best fitting model for initial propagation velocity
#=
Base model: $U = Γ/(4πR) * (log(8/α) - 1/2 + A)$
Expectation: fit A to some combination of $At$, $Re$, and $t$.
α=a/R is the aspect ratio of the vortex core.
Multiple choices viable for $a$, $R$, and $U$:
$a$:
 - per-core from $a^2= (R₂² - R²)
 - averaged from all core measurements for a given set of conditions
 - or estimated from $a = sqrt(4*ν*t)$
$R$:
 - per-core from vorticity moment integrals
 - taken from PIV core picking
 - or taken from PLIF core picking
$U$:
 - from $U₀exp(-kt)$ fit
 - or local linear fit

Evaluate all permutations of these choices and find the best fit.
=#
IClist = filter(r -> !r.shockrun, good_PIV)

IClist.U_fit = map(eachrow(IClist)) do runmeta
    zfits = Dict{Symbol, Vector{Float64}}() #Vector{Float64}[]
    cores = runmeta.cores
    for lr in (:left, :right)
        t = corefield(cores, lr, :preshock, :t)
        isempty(t) && continue
        x = corefield(cores, lr, :preshock, :x)
        z = corefield(cores, lr, :preshock, :z)

        @. z_c(i, p) = p[1]*(1 - exp(-p[2]*(i - p[3])))
        zfit = curve_fit(z_c, [0; t], [0; z], [0.5, 200, 1e-3])
        # push!(zfits, zfit.param)
        zfits[Symbol(lr)] = zfit.param
    end
    length(zfits) == 0 && return missing
    mean(zfit -> zfit[1]*zfit[2]*exp(-zfit[2]*runmeta.t_TSI), values(zfits))
end
##
gbIC = groupby(IClist, [:MST_gas, :MST_psig])
transform!(gbIC, x -> (; a_med = median([getfield.(x.∫ₗ, :a); getfield.(x.∫ᵣ, :a)])))
transform!(gbIC, x -> (; a_formula = sqrt(4*MST_state(x[1,:]).driven.nu*x[1,:].t_TSI)))

IClist.R_TSIpick = map(eachrow(IClist)) do runmeta
    zfits = Dict{Symbol, Vector{Float64}}() #Vector{Float64}[]
    core = load(datadir("PIV", "cores", runname(runmeta)*".jld2"))
    isnan(core["leftcore"][1]) || isnan(core["rightcore"][1]) && return nothing
    R = 1e-3*hypot((core["leftcore"] .- core["rightcore"])...) ./ 2
end

IClist.R_∫rω = map(eachrow(IClist)) do r
    (r.∫ᵣ.R - r.∫ₗ.R) / 2
end

##
IClist.U_predict_naive = map(eachrow(IClist)) do r
    R = r.R_TSIpick
    # R = ustrip(u"m", r.R_∫rω)
    α = ustrip(u"m", r.a_med) / R
    ν = MST_state(r).driven.nu
    t = r.t_TSI
    r.Γ / (4π*R*u"m") * (log(8/α) - 0.5)
end

IClist.U_predict = map(eachrow(IClist)) do r
    R = r.R_TSIpick
    # R = ustrip(u"m", r.R_∫rω)
    α = ustrip(u"m", r.a_med) / R
    ν = MST_state(r).driven.nu
    t = r.t_TSI
    r.Γ / (4π*R*u"m") * (log(8/α) - 0.5 - 1.6(0.6 - r.At)^2)
end
# Plot U_predict & U_fit versus time
f = Figure(resolution=100 .* (7, 6))

colorset = Dict(unique(IClist.MST_psig) .=> Makie.wong_colors()[1:7])
for i in 1:2, j in 1:2
    idx = j + 2*(i-1)
    gas, _ = shockICs[idx]
    local ax = Axis(f[i, j], xlabel="time [ms]", ylabel="u [m/s]",)
    for k in filter(k -> k.MST_gas == gas, sort(keys(gbIC), by=k -> k.MST_psig, rev=true))
        grp = gbIC[k]
        for r in eachrow(grp)
            U₀ = ustrip(u"m/s", MST_state(r).u2)
            scatterlines!(ax, 1000 .* [r.t_TSI, r.t_TSI], 
                [ustrip.(u"m/s", r.U_predict), r.U_fit], 
                marker=[:diamond, :x], color=colorset[k.MST_psig], label=@sprintf("%.1f", U₀))
        end
    end
    axislegend(gas*"\n uₚ [m/s]", merge=true, rowgap=1)
end
save(plotsdir("circulation_jump", "U_predict_fit_versus_time.svg"), f)
f

## Plot relative prediction error (U_predict - U_fit) / U_fit versus U_fit
f = Figure(resolution=100 .* (7, 6))

for i in 1:2, j in 1:2
    idx = j + 2*(i-1)
    gas, _ = shockICs[idx]
    local ax = Axis(f[i, j], xlabel="u [m/s]", ylabel="% relative error",)
    for k in filter(k -> k.MST_gas == gas, sort(keys(gbIC), by=k -> k.MST_psig, rev=true))
        grp = gbIC[k]
        for r in eachrow(grp)
            U₀ = ustrip(u"m/s", MST_state(r).u2)
            Uf = r.U_fit
            Upredict = ustrip.(u"m/s", r.U_predict)
            scatter!(ax, Uf, 100(Upredict - Uf) ./ Uf, 
                marker=:diamond, color=colorset[k.MST_psig], label=@sprintf("%.1f", U₀))
            hlines!(ax, [0], color=:black)
        end
    end
    xlims!(ax, nothing, ((i, j) == (1, 2) ? 120 : 100))
    axislegend(gas*"\n uₚ [m/s]", merge=true, rowgap=0, position=:rt)
end
save(plotsdir("circulation_jump", "predict_error_versus_U_fit.svg"), f)
f
## Predict pre-shock circulation for shocked runs
PSlist = filter(r -> r.shockrun, good_PIV)

PSlist.U_fit = map(eachrow(PSlist)) do runmeta
    zfits = runmeta.zfits
    length(zfits) == 0 && return missing
    mean(zfit -> zfit[1]*zfit[2]*exp(-zfit[2]*runmeta.t_TSI), values(zfits))
end
##
gbPS = groupby(PSlist, [:MST_gas, :MST_psig])
transform!(gbPS, x -> (; a_med = median([getfield.(x.∫ₗ, :a); getfield.(x.∫ᵣ, :a)])))

PSlist.R_TSIpick = map(eachrow(PSlist)) do runmeta
    core = load(datadir("PIV", "cores", runname(runmeta)*".jld2"))
    isnan(core["leftcore"][1]) || isnan(core["rightcore"][1]) && return nothing
    R = 1e-3*hypot((core["leftcore"] .- core["rightcore"])...) ./ 2
end

PSlist.Γ_pred = map(eachrow(PSlist)) do r
    R = r.R_TSIpick
    α = ustrip(u"m", r.a_med) / R
    r.U_fit*u"m/s" * 4π*R*u"m" / (log(8/α) - 0.5 - 1.6*(0.6 - r.At)^2)
end

## Plot circulation jump (relative to pre-shock prediction) versus time relative to t_SVI
f = Figure(resolution=100 .* (7, 5))
ax1 = Axis(f[1, 1], xlabel="post-shock time [ms]", ylabel="ΔΓ [m²/s]", 
        xminorticksvisible=true, yminorticksvisible=true)
ΔΓpred = Dict{String,Vector{Float64}}()
for i in 1:2, j in 1:2
    idx = j + 2*(i-1)
    gas, psig = shockICs[idx]
    shocks = gbPS[(MST_gas = gas, MST_psig=psig)]
    ΔΓ = shocks.Γ .- shocks.Γ_pred
    ΔΓpred[gas] = ustrip.(u"m^2/s", ΔΓ)
    plot!(ax1, 1e3*(shocks.t_TSI .- shocks.t_SVI), ustrip.(u"m^2/s", ΔΓ), 
        marker=:dtriangle, label=gas, color=Makie.wong_colors()[idx])
    hlines!(ax1, ustrip(u"m^2/s", mean(ΔΓ)), color=Makie.wong_colors()[idx], linestyle=:dot, label=gas)
end
axislegend(position = :rt, merge=true)
# Plot ΔΓ distribution versus Atwood number
ax2 = Axis(f[1, 2], xlabel="Atwood number", 
        xminorticksvisible=true, yminorticksvisible=true)
for i in 1:2, j in 1:2
    idx = j + 2*(i-1)
    gas, psig = shockICs[idx]
    ΔΓ = ΔΓpred[gas]
    shocks = gbPS[(MST_gas = gas, MST_psig=psig)]
    boxplot!(ax2, fill(mean(shocks.At), length(ΔΓ)), ΔΓ, 
        label=gas, color=Makie.wong_colors()[idx], width=0.07)
end
linkyaxes!(ax1, ax2)
save(plotsdir("circulation_jump", "Delta_Gamma_predicted_versus_At.svg"), f,)
f

## Plot normalized predicted and measured circulation jump versus Atwood number
f = Figure(resolution=100 .* (7, 5))
ax1 = Axis(f[1, 1], xlabel="Atwood number", ylabel=L"\Delta\Gamma/\Gamma_{-}", 
        xminorticksvisible=true, yminorticksvisible=true)


##
runmeta = select_run(runlist, "2023-01-19_run9")
f = Figure()
ax = Axis(f[1, 1])
scatter!(ax, getfield.(values(runmeta.cores[:left][:preshock]), :t),
                getfield.(values(runmeta.cores[:left][:preshock]), :z), color=:black, label="left")
scatter!(ax, getfield.(values(runmeta.cores[:right][:preshock]), :t),
                getfield.(values(runmeta.cores[:right][:preshock]), :z), color=:red, label="right")
scatter!(ax, getfield.(values(runmeta.cores[:left][:postshock]), :t),
                getfield.(values(runmeta.cores[:left][:postshock]), :z), color=:black, label="")
scatter!(ax, getfield.(values(runmeta.cores[:right][:postshock]), :t),
                getfield.(values(runmeta.cores[:right][:postshock]), :z), color=:red, label="")
vlines!(ax, [runmeta.t_SVI], color=:blue, label="SVI")
vlines!(ax, [runmeta.t_TSI], color=:black, label="TSI")
f