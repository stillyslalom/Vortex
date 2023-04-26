using DrWatson
@quickactivate "Vortex"
using GLMakie, JLD2, ImageFiltering, Unitful

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

##
savefigs = true
# begin; runmeta = select_run(runlist, "2022-12-20_run2")
foreach(eachrow(runlist)) do runmeta
    PIV = load(datadir("PIV", "registered", runname(runmeta)*".jld2"))
    ω = vorticity(PIV["x"] .* 1e-3, PIV["y"] .* 1e-3, PIV["u"], PIV["v"])
    core = load(datadir("PIV", "cores", runname(runmeta)*".jld2"))
    isnan(core["leftcore"][1]) || isnan(core["rightcore"][1]) && return nothing
    centerpoint = (core["leftcore"] .+ core["rightcore"]) ./ 2

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

    f = Figure()
    ax = Axis(f[1, 1], aspect=DataAspect(), xlabel="x [m]", ylabel="Z [m]", 
        title="Ring ROI for $(runname(runmeta))")
    hm = heatmap!(ax, ustrip.(u"m", xx .- xmid), ustrip.(u"m", yy), ustrip.(u"s^-1", ωᵥ),
        colormap=:balance, colorrange=(-5e4, 5e4))
    # Colorbar(f[1, 2], hm, label="ω [s⁻¹]", ticklabelsize=12)
    scatter!(ax, ustrip.(u"m", [∫ᵣ.R, ∫ₗ.R]), ustrip.(u"m", [∫ᵣ.Z, ∫ₗ.Z]), 
        markersize=40, marker='+', color=[:blue, :red])
    scatter!(ax, Circle(Point(ustrip.(u"m", ∫ᵣ.R), ustrip.(u"m", ∫ᵣ.Z)), ustrip.(u"m", ∫ᵣ.a)), 
        color=:blue, fill=false)
    scatter!(ax, Circle(Point(ustrip.(u"m", ∫ₗ.R), ustrip.(u"m", ∫ₗ.Z)), ustrip.(u"m", ∫ₗ.a)), 
        color=:red, fill=false)
    # resize_to_layout!(f)
    # rowsize!(f.layout, 1, Aspect(1, 1))
    savefigs && save(plotsdir("vortex_ROI", runname(runmeta)*".png"), f)
    f
end

##
# foreach(eachrow(runlist)) do runmeta
function vortex_integrals(runmeta)
    PIV = load(datadir("PIV", "registered", runname(runmeta)*".jld2"))
    ω = vorticity(PIV["x"] .* 1e-3, PIV["y"] .* 1e-3, PIV["u"], PIV["v"])
    core = load(datadir("PIV", "cores", runname(runmeta)*".jld2"))
    isnan(core["leftcore"][1]) || isnan(core["rightcore"][1]) && return nothing
    centerpoint = (core["leftcore"] .+ core["rightcore"]) ./ 2

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
    return (; leftcore = ∫ₗ, rightcore = ∫ᵣ)
end
# transform!(shockruns, AsTable(All()) => ByRow(PLIF_trajectory) => AsTable)

core_runlist = transform(runlist,  AsTable(All()) => ByRow(vortex_integrals) => AsTable)
##
core_runlist.shockrun = .!ismissing.(core_runlist.ptrace_path)
normalruns = filter(core_runlist) do row
    (20.5 .<= row.MST_psig .<= 25.5)
end

gdf = groupby(normalruns, [:shockrun, :MST_gas])

combine(gdf) do g
    Γₗ = median(getfield.(g.leftcore, :Γ))
    Γᵣ = median(getfield.(g.rightcore, :Γ))
    (; Γₗ, Γᵣ, Γ = (Γₗ - Γᵣ)/2)
end

gases = ["N2", "Ar", "CF4", "SF6"]
gasidx = Dict(gases .=> 1:4)

fig = Figure(resolution = (600, 400))
ax = Axis(fig[1, 1], xticks = (1:length(gases), gases), ylabel="Γ [m²/s]",
    title="Pre- and post-shock vortex circulation")
colors = categorical_colors(:Set1, length(gases))
for k in keys(gdf)
    Γₗ = getfield.(gdf[k].leftcore, :Γ)
    Γᵣ = getfield.(gdf[k].rightcore, :Γ)
    Γ = (Γₗ - Γᵣ)/2
    idx = gasidx[k.MST_gas]
    a = fill(idx - 0.15 + 0.3*k.shockrun, length(Γ))
    boxplot!(ax, a, ustrip.(Γ); whiskerwidth = 1, width = 0.35,
        color = (colors[idx], 0.35 + 0.4*k.shockrun), whiskercolor = (colors[idx], 1),
        mediancolor = :black)
end

fig