using DrWatson
@quickactivate "Vortex"
using Vortex
using TrixiShockTube
using DataFramesMeta
using LaTeXStrings
using Printf
using PrettyTables
using CSV
using Unitful
using PyThermo
using CairoMakie

runlist = loadmeta() do m
    !ismissing(m.TSI_ID) && !ismissing(m.registration_path) && !ismissing(m.cine_ID) && !ismissing(m.timings_path)
end

##
ICcondgrp = groupby(runlist, [:MST_gas, :MST_psig])

ICconds = combine(ICcondgrp, g -> (; IC = MST_state(first(eachrow(g)))))

ambient = Species("N2", P = 14.5u"psi", T=19u"°C")

shockICs = ["N2" => 20.5, "Ar" => 20.5, "CF4" => 25.5, "SF6" => 25.5,]

##
function run_MST(gas, psig)
    IC = MST_state(gas, psig)
    ambient = Species("N2", P = 14.5u"psi", T=19u"°C")
    slabs = (IC.driver => ustrip(u"m", 3.308u"inch"),
            IC.driven => ustrip(u"m", 10.61u"inch"), 
            ambient => ustrip(u"m", 30u"inch"))

    ic, equations = TrixiShockTube.build_shocktube_ic(slabs)
    semi = build_semidiscretization(slabs, ic, equations;
            initial_refinement_level=8)

    saveat = range(0, step=5e-6, stop=ustrip(u"s", 6slabs[2][2]*u"m" / (IC.Ms*soundspeed(IC.driven))))
    limiter! = build_limiter(equations)

    sol = run_shock(semi, saveat, 0.5, limiter!)

    x, t, data = xtdata(sol, semi)
end

transform!(ICconds, [:MST_gas, :MST_psig] => ByRow(run_MST) => [:x, :t, :data])

##
ICconds.u_p = map(eachrow(ICconds)) do r
    iₓ_exit = findfirst(>(0.3535172), r.x)::Int
    r.data[:v1][iₓ_exit, :]
end

ICconds.i0 = map(eachrow(ICconds)) do r
    i0 = findfirst(>(maximum(r.u_p)/10), r.u_p)::Int
end

ICconds.t0 = map(r -> r.t[r.i0], eachrow(ICconds))

ICconds.t = map(eachrow(ICconds)) do r
    r.t .- r.t0
end

ICconds.ifinal = map(eachrow(ICconds)) do r
    findfirst(<(0), r.u_p)::Int - 1
end

ICconds.tf = map(r -> r.t[r.ifinal], eachrow(ICconds))

##

ICconds.L_D = map(eachrow(ICconds)) do r
    cumsum(r.u_p .* (r.t[2] - r.t[1])) ./ ustrip(u"m", 0.875u"inch")
end

ICconds.i_pinchoff = map(eachrow(ICconds)) do r
    findlast(>(0.7*maximum(r.u_p)), r.u_p)::Int
end

ICconds.t_pinchoff = map(r -> r.t[r.i_pinchoff], eachrow(ICconds))

ICconds.L_D_pinchoff = map(eachrow(ICconds)) do r
    r.L_D[r.i_pinchoff]
end
##
f = Figure(resolution=100 .* (6, 5))
ax1 = Axis(f[1, 1], xticklabelsvisible=false, ylabel=L"$u_p$ [m/s]", 
    xminorticksvisible=true, yminorticksvisible=true)
ax2 = Axis(f[2, 1], ylabel=L"L/D", xlabel=L"$t$ [ms]", 
    xminorticksvisible=true, yminorticksvisible=true)
ICgrp = groupby(ICconds, [:MST_gas, :MST_psig])
gas_str = Dict("N2" => "N₂", "Ar" => "Ar", "CF4" => "CF₄", "SF6" => "SF₆")

for (gas, psig) in shockICs
    lines!(ax1, 1e3 .* only(ICgrp[(MST_gas = gas, MST_psig = psig)]).t, 
        only(ICgrp[(MST_gas = gas, MST_psig = psig)]).u_p,
        label = @sprintf("%s (%0.2f)", gas_str[gas], (psig + 14.5)/14.5))
    lines!(ax2, 1e3 .* only(ICgrp[(MST_gas = gas, MST_psig = psig)]).t, 
        only(ICgrp[(MST_gas = gas, MST_psig = psig)]).L_D)
    vlines!(ax1, 1e3 .* only(ICgrp[(MST_gas = gas, MST_psig = psig)]).t_pinchoff,
        linestyle=:dashdot)
    vlines!(ax2, 1e3 .* only(ICgrp[(MST_gas = gas, MST_psig = psig)]).t_pinchoff,
        linestyle=:dashdot)
end
xlims!(ax1, -1e-2, 2.25)
xlims!(ax2, -1e-2, 2.25)
ylims!(ax1, 0, 120)
ylims!(ax2, 0, 6)
f.current_axis[] = ax1
axislegend(L"Gas $(p_4/p_1)$", rowgap=1)
# save(plotsdir("discharge_profiles", "shockconds.svg"), f)
f
##

ICconds.p4_p1 = map(eachrow(ICconds)) do r
    pressure(r.IC.driver) / pressure(r.IC.driven)
end

ICconds.p4_kPa = map(eachrow(ICconds)) do r
    ustrip(u"kPa", pressure(r.IC.driver))
end

ICconds.Atwood = map(eachrow(ICconds)) do r
    ρ1 = PyThermo.density(r.IC.driven)
    ρ2 = PyThermo.density(Species("N2", T = temperature(r.IC.driven), P = pressure(r.IC.driven)))
    (ρ1 - ρ2) / (ρ1 + ρ2)
end

ICconds.ace_frac = map(eachrow(ICconds)) do r
    r.IC.driven.zs[2]
end

ICconds.Re_discharge = map(eachrow(ICconds)) do r
    mean(r.u_p[r.i0:r.i_pinchoff]) * ustrip(u"m", 0.875u"inch") / r.IC.driven.nu 
end

ICconds.driven_density = map(eachrow(ICconds)) do r
    PyThermo.density(r.IC.driven)
end

exportICs = select(ICconds, [:MST_gas, :MST_psig, :p4_kPa, :t0, :L_D_pinchoff, :p4_p1, :Atwood, :Re_discharge, :ace_frac, :driven_density])
exportICs.t0_μs = map(r -> r.t0 * 1e6, eachrow(exportICs))

sort!(exportICs, [:driven_density, :MST_psig])
# delete driven density column
select!(exportICs, Not(:driven_density))
select!(exportICs, Not(:t0))
CSV.write(plotsdir("discharge_profiles", "ICconds.csv"), exportICs)

open(plotsdir("discharge_profiles", "ICconds.tex"), "w") do io
    gas_strs = Dict("N2" => L"\mathrm{N_2}", "Ar" => L"\mathrm{Ar}", "CF4" => L"\mathrm{CF_4}", "SF6" => L"\mathrm{SF_6}")
    pretty_table(io, exportICs, backend=Val(:latex),
        header = ["MST gas", L"$p_4$ [psig]", L"$p_4$ [kPa]", L"$L/D$", L"p_4/p_1", L"A", L"Discharge $Re$", "% acetone", L"$\Delta t$ [\mu s]"],
            formatters = (v, i, j) -> begin
                if j == 1
                    return gas_strs[v]
                else
                    return round(v, sigdigits=2)
                end
            end
    )
end

## Make labeled x-t diagram for 20.5 psig Ar
ICconds_20p5 = only(groupby(ICconds, [:MST_gas, :MST_psig])[(MST_gas = "Ar", MST_psig = 20.5)])
i_x0 = findfirst(ICconds_20p5.data[:rho2] .> 1e-1)

n_pre = 150
p_pre = repeat(ICconds_20p5.data[:p][:,1], 1, n_pre)
p = hcat(p_pre, ICconds_20p5.data[:p])  

x0 = ICconds_20p5.x[i_x0]
xx = ICconds_20p5.x .- x0

Δt = mean(diff(ICconds_20p5.t))
tt = range(stop = ICconds_20p5.t[end], length = size(p, 2), step = Δt)

f = Figure(resolution=(600, 400))
ax = Axis(f[1, 1], xlabel = L"$x$ [m]", ylabel = L"$t$ [ms]", 
    xminorticksvisible=true, yminorticksvisible=true)
hm = heatmap!(ax, xx, 1e3 .* tt, p ./ 1e5, 
    colormap = :viridis, colorrange = (1, 2.5))
contour!(ax, xx, 1e3 .* ICconds_20p5.t, ICconds_20p5.data[:rho2], 
    levels = [1e-1], 
    linewidth = 1.5, color = :black, linealpha = 0.5)
xlims!(ax, nothing, 0.45)
ylims!(ax, -0.18, 1.2)
Colorbar(f[1, 2], hm, label = L"$p$ [bar]")
text!(0.01 - x0, 0.15; text="🅡", fontsize=24)
text!(0.015 - x0, 0.04; text="❹", fontsize=24)
text!(0.06 - x0, 0.13; text="❸", fontsize=24, color=:white)
text!(0.105 - x0, 0.13; text="❷", fontsize=24, color=:white)
text!(0.15 - x0, 0.04; text="❶", fontsize=24, color=:white)

bracket!(ax, xx[1], 0, 0, 0, text="Driver", orientation=:down)
x_ambient = xx[findfirst(>(0.1), ICconds_20p5.data[:rho3])]
bracket!(ax, 0, 0, x_ambient, 0, text="Driven", orientation=:down, 
    color=:white, textcolor=:white)
bracket!(ax, x_ambient, 0, 0.45, 0, text="Ambient", orientation=:down, 
    color=:white, textcolor=:white)
save(plotsdir("discharge_profiles", "xt_diagram.png"), f)
f
