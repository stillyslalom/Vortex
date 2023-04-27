using DrWatson
@quickactivate "Vortex"
using Vortex
using JLD2
using HDF5
using PyThermo
using Unitful
using Roots
using OrderedCollections, FileIO
using Dates
using Tar, CodecZlib

##
# HDF5 schema:
# grid_x::Vector{Float64}
# grid_y::Vector{Float64}
# time::Vector{Float64}
# /PLIF
#    001::Matrix{Float32}
#    002::Matrix{Float32}
#    ...
# /PIV
#    u::Matrix{Float64}
#    v::Matrix{Float64}
#    sigma_u::Matrix{Float64}
#    sigma_v::Matrix{Float64}
#    status::Matrix{Int}
# /metadata
#    ring_gas::String
#    discharge_velocity::Float64
#    i_PIV::Int
#    run_ID::String

runlist = loadmeta(m -> !ismissing(m.TSI_ID) && !ismissing(m.registration_path) 
                     && !ismissing(m.cine_ID) && !ismissing(m.timings_path))
Vortex.add_registrations!(runlist)

qualities = map(eachrow(runlist)) do runmeta
    qualities = JLD2.load(datadir("PIV", "masks", runname(runmeta) * ".jld2"))
    qualities["runname"] = runname(runmeta)
    qualities
end |> DataFrame

translations = map(eachrow(runlist)) do runmeta
    translations = JLD2.load(datadir("PIV", "manual_translations", runname(runmeta) * ".jld2"))
    translations["runname"] = runname(runmeta)
    translations
end |> DataFrame

runlist.runname = runname.(eachrow(runlist))
leftjoin!(runlist, DataFrame(qualities), on=:runname)
leftjoin!(runlist, DataFrame(translations), on=:runname)
filter!(row -> row.quality == 2 && ismissing(row.ptrace_path), runlist)
filter!(r -> !(Date(2022,12,15) <= r.Date <= Date(2022, 12, 16)), runlist)

##
function MST_state(runmeta)
    p_ace = Species("acetone", T = 22u"°C").Psat*u"Pa"
    p_total = (runmeta.MST_psig + 14.5)*u"psi" |> u"Pa"
    χ_ace = p_ace/p_total
    MSTgas = Mixture([runmeta.MST_gas => 1-χ_ace, "acetone" => χ_ace])
    Ma = fzero(1.31) do M
        pressure(PyThermo.ShockTube.driverpressure(MSTgas, MSTgas, M)) - p_total
    end
    sc = PyThermo.ShockTube.shockcalc(MSTgas, MSTgas, Ma)
end

##
# begin; runmeta = rand(eachrow(runlist))
foreach(eachrow(runlist)) do runmeta
    cine = Vortex.phantom_bgsub(runmeta)
    i_PIV, t_err = Vortex.phantom_timing(runmeta, Vortex.read_timings(runmeta), size(cine, 3))
    t_err > 0 && @warn "PIV time is not simultaneous with cine frame"
    cine_PIV = cine[:, :, i_PIV]
    uv = Vortex.transform_prana(runmeta, axes(cine_PIV); tx=-runmeta.tx, ty=runmeta.ty)
    phantomscale = runmeta.phantomscale
    xx = range(0, step=phantomscale.mm_per_px, length=size(cine_PIV, 2))
    yy = range(stop = phantomscale.origin_mm, step=phantomscale.mm_per_px, length=size(cine_PIV, 1))
    u, v = (reverse.((uv.u, uv.v), dims=2)...,)
    σu, σv, status = (reverse.((uv.σu, uv.σv, uv.status), dims=2)...,)

    sc = MST_state(runmeta)
    Ws = soundspeed(sc.driven)*sc.Ms # shock speed in driven gas
    t_shockprop = 10.5u"inch" / Ws

    timings = read_timings(runmeta)
    rupture_lag = (timings[:Rupture1].Shock + timings[:Rupture1].Delay)*u"μs" + t_shockprop
    # First Phantom frame occurs at t0_Phantom after rupture
    t0_Phantom = timings[:Death_Star].Delay - ustrip(u"μs", rupture_lag)

    t_Phantom = range(t0_Phantom + 50, length=size(cine, 3), step=50) ./ 1000
    data = OrderedDict{String, Any}(
        "grid_x" => collect(xx ./ 1000), # mm to m
        "grid_y" => collect(yy ./ 1000), # mm to m
        "time" => collect(t_Phantom ./ 1e6), # μs to s
    )
    foreach(1:size(cine, 3)) do i
        data["PLIF/"*lpad(i, 3, '0')] = rotr90(@view cine[:, :, i])
    end
    data["PIV/u"] = u
    data["PIV/v"] = v
    data["PIV/sigma_u"] = σu
    data["PIV/sigma_v"] = σv
    data["PIV/status"] = status
    data["metadata/ring_gas"] = runmeta.MST_gas
    data["metadata/discharge_velocity"] = ustrip(u"m/s", sc.u2)
    data["metadata/i_PIV"] = i_PIV
    data["metadata/run_ID"] = runname(runmeta)

    save(datadir("export", runname(runmeta) * ".h5"), data)
    nothing
end

##
open(datadir("export.tar.gz"), "w") do tar_gz
    tar = GzipCompressorStream(tar_gz)
    Tar.create(datadir("export"), tar)
end
