using DrWatson
@quickactivate "Vortex"
using Vortex, StaticArrays, CoordinateTransformations, Interpolations, ImageTransformations
using GLMakie
using Interpolations: Line

# Load runs with both PLIF & PIV data
runlist = loadmeta() do meta
    !ismissing(meta.TSI_ID) && !ismissing(meta.cine_ID) && !ismissing(meta.timings_path)
end

Vortex.add_registrations!(runlist)

##
# Load & transform data for a random run
runmeta = rand(eachrow(runlist))
pranaraw = PranaData(datadir("PIV", "runs", runname(runmeta)))
PIV_pixgrid = (pranaraw[end].x, pranaraw[end].y) ./ runmeta.dx
corners = [SVector(PIV_pixgrid[1][1], PIV_pixgrid[2][1]),
           SVector(PIV_pixgrid[1][1], PIV_pixgrid[2][end]),
           SVector(PIV_pixgrid[1][end], PIV_pixgrid[2][end]),
           SVector(PIV_pixgrid[1][end], PIV_pixgrid[2][1])]

M = LinearMap(SMatrix{3,3,Float32}(runmeta.transform'))
tf = Vortex.TSI2Phantom(M)
invtf = Vortex.Phantom2TSI(M)
PIV_corners = map(reverse ∘ tf ∘ reverse, corners)
# bounding box of PIV data in Phantom pixel coordinates
jmin, jmax = extrema(getindex.(PIV_corners, 1))
imin, imax = extrema(getindex.(PIV_corners, 2))
irange = range(round.(Int, (imin, imax))...)
jrange = range(round.(Int, (jmin, jmax))...)

# u_itp = linear_interpolation(PIV_pixgrid, pranaraw[end].u, extrapolation_bc=Line())
# v_itp = linear_interpolation(PIV_pixgrid, pranaraw[end].v, extrapolation_bc=Line())
u_itp = linear_interpolation(reverse(PIV_pixgrid), 
                             collect(reverse(pranaraw[end].u, dims=1)'),
                             extrapolation_bc=Line())
v_itp = linear_interpolation(reverse(PIV_pixgrid),
                             collect(reverse(pranaraw[end].v, dims=1)'),
                             extrapolation_bc=Line())

u_Phantom = warp(u_itp, invtf, (irange, jrange))
v_Phantom = warp(v_itp, invtf, (irange, jrange))
heatmap(parent(v_Phantom), axis=(; aspect=DataAspect()))

##
timings = read_timings(runmeta)
t_TSI = (timings[:PIV_trig].Delay + runmeta.PIV_delay) / 1e6
t0_Phantom = timings[:Death_Star].Delay + 50
cine = phantom_bgsub(runmeta)
t_Phantom = range(t0_Phantom, length=size(cine, 3), step=50) / 1e6
t_err, cine_PIV_idx = findmin(t -> abs(t - t_TSI), t_Phantom)
cine_PIV = cine[:,:,cine_PIV_idx]
# f, ax, hm = heatmap(imadjust(rotr90(cine_PIV)), colormap=:grays, axis=(; aspect=DataAspect(),))
f, ax, hm = image(collect(cine_PIV'), axis=(aspect=DataAspect(), yreversed=true))
U_Phantom_itp = linear_interpolation((irange, jrange), parent(u_Phantom), extrapolation_bc=Line())
V_Phantom_itp = linear_interpolation((irange, jrange), parent(v_Phantom), extrapolation_bc=Line())
V̄ = _ -> 0 #linear_interpolation(irange, dropdims(median(parent(v_Phantom), dims=2), dims=2), extrapolation_bc=Line())
u(x, y) = Point2(V_Phantom_itp(x, y), -U_Phantom_itp(x, y) - V̄(y))
streamplot!(ax, u, jrange, irange, stepsize=1, colormap=:RdBu, linewidth=0.5)
# heatmap!(ax, jrange, irange, parent(u_Phantom), alpha=0.1)
scatter!(ax, PIV_corners, color=:red)
f

##
cine2prana(M, cine) = PerspectiveMap() ∘ inv(M) ∘ push1 ∘ AffineMap(SA[0 1; -1 0], SA[0, size(cine, 1) + 1])
itf = cine2prana(M, cine_PIV)
u_itp = linear_interpolation(PIV_pixgrid, pranaraw[end].u, extrapolation_bc=NaN)
heatmap(warp(u_itp, itf, (1:size(cine_PIV, 1), 1:size(cine_PIV, 2))) |> parent; axis=(; aspect=DataAspect()))

##
savefigs=true
# foreach(eachrow(runlist)) do runmeta
let runmeta = select_run(runlist, "2023-01-23_run1")

    pranaraw = PranaData(datadir("PIV", "runs", runname(runmeta)))
    GOOD = Vortex.MedianFilter(3, 10, Vortex.MedianComponents())(pranaraw[3])
    u′ = copy(pranaraw[end].u)
    v′ = copy(pranaraw[end].v)
    @. u′[!GOOD] = NaN
    @. v′[!GOOD] = NaN

    M = LinearMap(SMatrix{3,3,Float32}(runmeta.transform'))

    timings = read_timings(runmeta)
    t_TSI = (timings[:PIV_trig].Delay + runmeta.PIV_delay) / 1e6
    t0_Phantom = timings[:Death_Star].Delay + 50
    cine = phantom_bgsub(runmeta)
    t_Phantom = range(t0_Phantom, length=size(cine, 3), step=50) / 1e6
    t_err, cine_PIV_idx = findmin(t -> abs(t - t_TSI), t_Phantom)
    cine_PIV = cine[:,:,cine_PIV_idx]

    tf = Vortex.TSI2Phantom(M)
    itf = Vortex.Phantom2TSI(M)

    # Find warped dimensions of TSI pixels in Phantom space
    PLIFx = range(0, length=size(cine_PIV, 2), step=runmeta.phantomscale.mm_per_px)
    PLIFz = range(start=runmeta.phantomscale.origin_mm, length=size(cine_PIV, 1), step=-runmeta.phantomscale.mm_per_px)
    PLIFx_itp = linear_interpolation(axes(PLIFx), PLIFx, extrapolation_bc=Line())
    PLIFz_itp = linear_interpolation(axes(PLIFz), PLIFz, extrapolation_bc=Line())
    PLIFxz(i, j) = Point2(PLIFx_itp(j), PLIFz_itp(i))
    PLIFxz(p::Point2) = Point2(PLIFx_itp(p[2]), PLIFz_itp(p[1]))
    jj, ii = (pranaraw[end].x, pranaraw[end].y) ./ runmeta.dx
    Δxz(i, j, δi, δj) = PLIFxz(tf(Point2(i + δi, j + δj))) - PLIFxz(tf(Point2(i, j)))
    u = Array{Float64}(undef, length(ii), length(jj))
    v = Array{Float64}(undef, length(ii), length(jj))
    u_scale = [Δxz(i, j, 0, 1) for i in ii, j in jj]
    v_scale = [Δxz(i, j, -1, 0) for i in ii, j in jj]
    u = map(CartesianIndices(u′)) do idx
        i_prana, j_prana = idx.I
        j = i_prana
        i = size(u′, 2) - j_prana + 1
        (u′[i_prana, j_prana] * u_scale[i, j][1] + v′[i_prana, j_prana] * v_scale[i, j][1]) * 1000 / runmeta.dx
    end
    v = map(CartesianIndices(v′)) do idx
        i_prana, j_prana = idx.I
        j = i_prana
        i = size(v′, 2) - j_prana + 1
        (u′[i_prana, j_prana] * u_scale[i, j][2] + v′[i_prana, j_prana] * v_scale[i, j][2]) * 1000 / runmeta.dx
    end

    u_itp = linear_interpolation(reverse(PIV_pixgrid), reverse(u, dims=2)', extrapolation_bc=NaN)
    v_itp = linear_interpolation(reverse(PIV_pixgrid), reverse(v, dims=2)', extrapolation_bc=NaN)
    uw = collect(warp(u_itp, itf, axes(cine_PIV))')
    vw = collect(warp(v_itp, itf, axes(cine_PIV))')
    U = hypot.(uw, vw)
    Ū = map(eachcol(U)) do u
        u_filt = filter(!isnan, u)
        length(u_filt) == 0 ? NaN : median(u_filt)
    end
    V̄ = median(filter(!isnan, vw))
    Ū_itp = linear_interpolation(axes(Ū), Ū)
    u_phantom = linear_interpolation(axes(uw), uw, extrapolation_bc=NaN)
    v_phantom = linear_interpolation(axes(vw), vw, extrapolation_bc=NaN)
    u(x, y) = Point2(u_phantom(x, y), -v_phantom(x, y) + V̄)# - Ū_itp(y))
    f = Figure(resolution=(800, 1000))
    ax = Axis(f[1, 1]; aspect=DataAspect(), yreversed=true)
    hm = image!(ax, imadjust(cine_PIV', qmax=0.9995))
    # heatmap!(ax, collect(warp(v_itp, itf, axes(cine_PIV))'), colormap=[RGBA(0,0,1,0), RGBA(1,0,0,1)])
    # heatmap!(ax, U, colormap=[RGBA(0,0,0,0), RGBA(1,0,0,0.3)], colorrange=quantile(filter(!isnan, U), (0.01, 0.99)))
    streamplot!(ax, u, jrange, irange, stepsize=5, colormap=[RGBA(0,0.2,1,1), RGBA(1, 1, 0, 1)], linewidth=1,
    colorrange=quantile(filter(!isnan, U), (0.01, 0.99)), arrow_size=5, gridsize=(64, 64, 64))
    # arrows!(ax, axes(uw)..., uw, -vw, color=:red, lengthscale=0.01, arrowsize=3)

    # hidexdecorations!(ax)
    # hideydecorations!(ax)
    # ax2 = Axis(f[1, 1], aspect=DataAspect(), xlabel="x (mm)", ylabel="z (mm)",
    #     title = string(runname(runmeta), " - ", runmeta.MST_gas, " at ", runmeta.MST_psig, "psig"))
    # @unpack origin_mm, mm_per_px = runmeta.phantomscale
    # z_Phantom = range(stop=origin_mm, step=mm_per_px, length=size(U, 2))
    # x_Phantom = range(0, step=mm_per_px, length=size(U, 1))
    # heatmap!(ax2, x_Phantom, z_Phantom, fill(NaN, size(U)))

    savefigs && save(plotsdir("PLIF_PIV", runname(runmeta) * ".png"), f)
    f
end

## PIV skew correction scratch space
runmeta = select_run(runlist, "2023-01-19_Ar_IC1")
pranaraw = PranaData(datadir("PIV", "runs", runname(runmeta)))
M = LinearMap(SMatrix{3,3,Float32}(runmeta.transform'))
tf = Vortex.TSI2Phantom(M)
itf = Vortex.Phantom2TSI(M)
PLIFx = range(0, length=size(cine_PIV, 2), step=runmeta.phantomscale.mm_per_px)
PLIFz = range(start=runmeta.phantomscale.origin_mm, length=size(cine_PIV, 1), step=-runmeta.phantomscale.mm_per_px)
PLIFx_itp = linear_interpolation(axes(PLIFx), PLIFx, extrapolation_bc=Line())
PLIFz_itp = linear_interpolation(axes(PLIFz), PLIFz, extrapolation_bc=Line())
PLIFxz(i, j) = Point2(PLIFx_itp(j), PLIFz_itp(i))
PLIFxz(p::Point2) = Point2(PLIFx_itp(p[2]), PLIFz_itp(p[1]))
jj, ii = (pranaraw[end].x, pranaraw[end].y) ./ runmeta.dx
Δxz(i, j, δi, δj) = PLIFxz(tf(Point2(i + δi, j + δj))) - PLIFxz(tf(Point2(i, j)))
u_scale = [Δxz(i, j, 0, 1) for i in ii, j in jj]
v_scale = [Δxz(i, j, -1, 0) for i in ii, j in jj]