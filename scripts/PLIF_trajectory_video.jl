using DrWatson
@quickactivate "Vortex"
using Vortex
using GLMakie
using JLD2
using LaTeXStrings
using Printf
using Unitful
using VideoIO
using Interpolations

update_theme!() #Theme(fonts = (; regular = "New Roman", bold = "New Roman Bold")))
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
# Select the Ar run with the most registered preshock cores
Ar_runlist = filter(r -> r.MST_gas == "Ar", runlist)
i_max, _ = findmax(eachrow(Ar_runlist)) do r
    cores = r.cores
    length(cores[:left][:preshock]) + length(cores[:right][:preshock])
end

runmeta = Ar_runlist[i_max, :]

## Create x(i) and z(i) interpolators for left and right cores
cores = runmeta.cores
function core_itp(cores, xz)
    ii = collect(keys(cores))
    vv = [getfield(cores[i], xz) for i in ii]
    LinearInterpolation(ii, vv, extrapolation_bc = Line())
end

lx = core_itp(cores[:left][:preshock], :x)
lz = core_itp(cores[:left][:preshock], :z)
rx = core_itp(cores[:right][:preshock], :x)
rz = core_itp(cores[:right][:preshock], :z)

## Load and plot a frame from the cine
cine = phantom_bgsub(runmeta)
xx, zz = runmeta.grid.x, runmeta.grid.z
tt = range(runmeta.t₀, step=50e-6, length=size(cine, 3))

f = Figure(resolution=(430,550))

i = Observable(80)
frame = @lift imadjust(rotr90(cine[:,:,$i]), qmin=0.4, qmax=0.9995)
timelabel = @lift @sprintf("t = %.2f ms", tt[$i]*1e3)
lcore = @lift Point2f(lx($i), lz($i))
rcore = @lift Point2f(rx($i), rz($i))
ltraj = Observable(Point2f[])
rtraj = Observable(Point2f[])
wtraj = @lift reverse([max(1, 3exp(-i/20)) for i in 1:length($ltraj)])

on(lcore) do c
    push!(ltraj[], c - Point2f(0, 0.003))
    notify(ltraj)
end
on(rcore) do c
    push!(rtraj[], c - Point2f(0, 0.003))
    notify(rtraj)
end

ax = Axis(f[1, 1], aspect=DataAspect(), xminorticksvisible=true, yminorticksvisible=true,
    xlabel = L"$x$ [m]", ylabel = L"$z$ [m]")
p = image!(ax, xx, zz, frame)
text!(ax, xx[1] + 0.01, zz[1] + 0.01, text=timelabel, 
    color=:white, font = "Consolas")

scatter!(ax, lcore, color=:blue, marker='○', markersize=20)
scatter!(ax, rcore, color=:red, marker='○', markersize=20)
lines!(ax, ltraj, color=:blue, linewidth=wtraj)
lines!(ax, rtraj, color=:red, linewidth=wtraj)

ylims!(ax, extrema(zz))

f
## Record animation
record(f, plotsdir("trajectory_PLIF_video", "traj.mp4"), 9:100;
        framerate = 15) do idx
    i[] = idx
end
