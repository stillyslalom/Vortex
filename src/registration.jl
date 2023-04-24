"Distance between two points"
dist2d(x, y) = sqrt(sum((x .- y).^2))
dist2d(x) = dist2d(x...)

function add_registrations!(runlist)
    calibrations = collect_results(datadir("calibrations"), verbose=false)
    calibrations.path = basename.(calibrations.path)
    leftjoin!(runlist, calibrations, on=:registration_path => :path, matchmissing=:equal)
end

push1(x) = CoordinateTransformations.push(x, 1)
TSI2Phantom(M) = PerspectiveMap() ∘ M ∘ push1
Phantom2TSI(M) = PerspectiveMap() ∘ inv(M) ∘ push1

PLIFxz(itps, i, j) = Point2(itps.x(j), itps.z(i))
PLIFxz(itps, p::Point2) = Point2(itps.x(p[2]), itps.z(p[1]))
Δxz(itps, tf, i, j, δi, δj) = PLIFxz(itps, tf(Point2(i + δi, j + δj))) - PLIFxz(itps, tf(Point2(i, j)))

function calculate_dewarp_scale(itps, tf, ii, jj)
    u_scale = [Δxz(itps, tf, i, j, 0, 1) for i in ii, j in jj]
    v_scale = [Δxz(itps, tf, i, j, -1, 0) for i in ii, j in jj]
    u_scale, v_scale
end

function dewarp_displacements(u′, v′, u_scale, v_scale, dx)
    uv = map(CartesianIndices(u′)) do idx
        i_prana, j_prana = idx.I
        j = i_prana
        i = size(u′, 2) - j_prana + 1
        u = (u′[i_prana, j_prana] * u_scale[i, j][1] + v′[i_prana, j_prana] * v_scale[i, j][1]) * 1000 / dx #runmeta.dx
        v = (u′[i_prana, j_prana] * u_scale[i, j][2] + v′[i_prana, j_prana] * v_scale[i, j][2]) * 1000 / dx #runmeta.dx
        (; u, v)
    end
    StructArray(uv)
end

"""
    transform_prana(runmeta, cine_axes)

Transforms the PIV displacements from the TSI coordinate system to the Phantom coordinate system.
"""
function transform_prana(runmeta, cine_axes)
    pranaraw = PranaData(datadir("PIV", "runs", runname(runmeta)))
    u′, v′, status = Vortex.load_infilled_PIV(runmeta)
    M = LinearMap(SMatrix{3,3,Float32}(runmeta.transform'))
    tf = Vortex.TSI2Phantom(M)
    itf = Vortex.Phantom2TSI(M)
    phantomscale = runmeta.phantomscale

    # Find warped dimensions of TSI pixels in Phantom space
    PLIFx = range(0, length=length(cine_axes[2]), step=runmeta.phantomscale.mm_per_px)
    PLIFz = range(start=phantomscale.origin_mm, length=length(cine_axes[1]), step=-phantomscale.mm_per_px)
    PLIFx_itp = linear_interpolation(axes(PLIFx), PLIFx, extrapolation_bc=Line())
    PLIFz_itp = linear_interpolation(axes(PLIFz), PLIFz, extrapolation_bc=Line())
    itps = (; x=PLIFx_itp, z=PLIFz_itp)
    jj, ii = (pranaraw[end].x, pranaraw[end].y) ./ runmeta.dx
    u_scale, v_scale = calculate_dewarp_scale(itps, tf, ii, jj)
    uv = dewarp_displacements(u′, v′, u_scale, v_scale, runmeta.dx)

    PIV_pixgrid = (pranaraw[end].x, pranaraw[end].y) ./ runmeta.dx
    u_itp = linear_interpolation(reverse(PIV_pixgrid), reverse(uv.u, dims=2)', extrapolation_bc=NaN)
    v_itp = linear_interpolation(reverse(PIV_pixgrid), reverse(uv.v, dims=2)', extrapolation_bc=NaN)
    uw = collect(warp(u_itp, itf, cine_axes)')
    vw = collect(warp(v_itp, itf, cine_axes)')
    return StructArray(u = uw, v = vw)
end
