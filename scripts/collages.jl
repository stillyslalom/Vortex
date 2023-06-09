using DrWatson
@quickactivate "Vortex"
using Vortex
using CairoMakie
using JLD2
using LaTeXStrings
using Printf
using Unitful
using VideoIO

update_theme!(Theme(fonts = (; regular = "New Roman", bold = "New Roman Bold")))
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

runlist.PIV = [load(datadir("PIV", "registered", r.runname*".jld2")) for r in eachrow(runlist)]
runlist.PIVcore = [load(datadir("PIV", "cores", r.runname*".jld2")) for r in eachrow(runlist)]
runlist.PIV_quality = map(eachrow(runlist)) do runmeta
    maskpath = datadir("PIV", "masks", runname(runmeta)*".jld2")
    !isfile(maskpath) && return missing
    load(maskpath, "quality")
end
##
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
f = Figure(resolution=(800,540))

runs = select_run.(Ref(runlist), ["2023-01-23_run4" "2023-01-20_run5"; "2023-01-11_run5" "2023-01-18_run6"])
pretty_gas = Dict("N2" => "N₂", "Ar" => "Ar", "CF4" => "CF₄", "SF6" => "SF₆")
axs = Dict()
for i in 1:2, j in 1:2
    local ax = Axis(f[i,j], aspect=DataAspect(), xminorticksvisible=true, yminorticksvisible=true,
        xlabel = L"x / R", ylabel=L"(z - z_c) / R", title=pretty_gas[runs[i,j].MST_gas])
    push!(axs, (i,j) => ax)
    core = runs[i,j].PIVcore
    PIV = runs[i,j].PIV
    xx = PIV["x"]
    zz = PIV["y"]
    ω = Vortex.vorticity(PIV["u"], PIV["v"], 1e-3*step(xx), 1e-3*step(zz))
    ∫ₗ = runs[i,j].∫ₗ
    ∫ᵣ = runs[i,j].∫ᵣ 

    R = mean((ustrip.(u"mm", [∫ᵣ.R, - ∫ₗ.R])))
    xc = (core["leftcore"][1] + core["rightcore"][1]) / 2
    zc = mean(ustrip.(u"mm", [∫ᵣ.Z, ∫ₗ.Z]))
    hm = heatmap!(ax, (xx .- xc) ./ R, (zz .- zc) ./ R, ω,
        colormap=:balance, colorrange = (-1e5, 1e5), rasterize=true)
    scatter!(ax, (ustrip.(u"mm", [∫ᵣ.R, ∫ₗ.R])) ./ R, (ustrip.(u"mm", [∫ᵣ.Z, ∫ₗ.Z]) .- zc) ./ R, 
        markersize=30, marker='+', color=[:aqua, :red],)
    lines!(ax, Circle(Point(ustrip.(u"mm", ∫ᵣ.R) ./ R, (ustrip.(u"mm", ∫ᵣ.Z) .- zc) ./ R), ustrip.(u"mm", ∫ᵣ.a) ./ R), 
        color=:aqua, fill=false, linewidth=3)
    lines!(ax, Circle(Point(ustrip.(u"mm", ∫ₗ.R) ./ R, (ustrip.(u"mm", ∫ₗ.Z) .- zc) ./ R), ustrip.(u"mm", ∫ₗ.a) ./ R), 
        color=:red, fill=false, linewidth=3)
    xlims!(ax, -1.5, 1.5)
    ylims!(ax, -1, 1)
end
hidexdecorations!.((axs[(1,1)], axs[(1, 2)]), ticks=false, minorticks=false)
hideydecorations!.((axs[(1, 2)], axs[(2, 2)]), ticks=false, minorticks=false)
Colorbar(f[1:2, 3], colormap=:balance, colorrange = (-1, 1), label=L"\omega / 10^5 \, [s^{-1}]")

save(plotsdir("collages", "vorticity.svg"), f)
f

## Montage of one frame from each PLIF video
savevids = false
if savevids
    cineframes = map(enumerate(eachrow(runlist))) do (i, m)
        phantom_bgsub(m)[:,:,m.i_TSI]
    end
    ##
    cineframes_adj = map(cineframes) do cf
        imadjust(cf; qmax=0.9995) .|> N0f8
    end
    ##  
    VideoIO.save(plotsdir("collages", "phantom_bgsub.mp4"), cineframes_adj, framerate=100)
    VideoIO.save(plotsdir("collages", "phantom_bgsub_slow.mp4"), cineframes_adj, framerate=20)
end

## Collage of each set of PIV vorticity fields sorted by post-shock time
testconds = Dict("N2" => 20.5, "Ar" => 20.5, "CF4" => 25.5, "SF6" => 25.5)
shockruns = filter(runlist) do r
    !ismissing(r.ptrace_path) && r.MST_gas ∈ keys(testconds) && r.MST_psig == testconds[r.MST_gas]
end

shockruns.t_rel_TSI = shockruns.t_TSI - shockruns.t_SVI
preshocks = filter(r -> r.t_rel_TSI < 0, shockruns)
filter!(r -> r.t_rel_TSI > 0 && r.PIV_quality == 2, shockruns)

sort!(shockruns, :t_rel_TSI)
shockgroups = groupby(shockruns, :MST_gas)

