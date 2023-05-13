using DrWatson
@quickactivate "Vortex"
using Vortex
using TrixiShockTube
using DataFramesMeta
using LaTeXStrings
using Printf

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
    iₓ_exit = findfirst(>(sum(last.(slabs)[1:2])), x)::Int
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
save(plotsdir("discharge_profiles", "shockconds.svg"), f)
f
##
