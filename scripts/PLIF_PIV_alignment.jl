using DrWatson
@quickactivate "Vortex"
using Colors, GLMakie

# Load runs with both PLIF & PIV data
runlist = loadmeta() do meta
    !ismissing(meta.TSI_ID) && !ismissing(meta.cine_ID) && !ismissing(meta.timings_path)
end

Vortex.add_registrations!(runlist)

##
function manual_translation(runmeta)
    cine = Vortex.phantom_bgsub(runmeta)
    i_PIV, t_err = Vortex.phantom_timing(runmeta, Vortex.read_timings(runmeta), size(cine, 3))
    t_err > 0 && @warn "PIV time is not simultaneous with cine frame"
    cine_PIV = cine[:, :, i_PIV]

    tx = Observable(0)
    ty = Observable(0)
    uv = Vortex.transform_prana(runmeta, axes(cine_PIV))

    f = Figure(resolution=(1100,1300))
    sg = SliderGrid(f[2, 1],
        (label = "X shift", range = -50:50, startvalue = 0),
        (label = "Y shift", range = -50:50, startvalue = 0),
        (label = "-log(1 - qmax)", range=0.5:0.1:10, startvalue=4))
    on(sg.sliders[1].value) do v
        tx[] = v[]
    end
    on(sg.sliders[2].value) do v
        ty[] = v[]
    end

    ax = Axis(f[1, 1], aspect=DataAspect(), yreversed=true)
    cine_adjusted = @lift imadjust(Float32.(cine_PIV), qmax=(1 - exp10(-$(sg.sliders[3].value))))'
    cine_ax1 = @lift axes($cine_adjusted, 1) .+ $tx
    cine_ax2 = @lift axes($cine_adjusted, 2) .- $ty
    heatmap!(ax, cine_ax1, cine_ax2, cine_adjusted, 
        colormap=[RGBA(1, 1, 1, 1), RGBA(0, 0, 0, 1)])

    v_shocked = filter(<(-200), uv.v)
    V₀ = isempty(v_shocked) ? 0 : median(v_shocked)
    v′ = copy(uv.v)
    v′[v′ .< V₀/2] .-= V₀
    function u_xy(u, v)
        u_itp = linear_interpolation(axes(u), u, extrapolation_bc=NaN)
        v_itp = linear_interpolation(axes(v), v, extrapolation_bc=NaN)
        (x, y) -> Point2(u_itp(x, y), -v_itp(x, y))
    end
    u = u_xy(uv.u, v′)

    streamplot!(ax, u, axes(v′)..., stepsize=0.5, gridsize=(128, 128, 128), arrow_size=1)

    wait(GLMakie.Screen(f.scene)) # wait for plot to close
    Dict("tx"=> tx[], "ty"=>ty[])
end

##
# runmeta = rand(eachrow(runlist))
translations = map(eachrow(runlist)) do runmeta
    data, path = produce_or_load(manual_translation, runmeta, datadir("PIV", "manual_translations");
                     filename=runname, tag=false)
    data["runname"] = runname(runmeta)
    data
end

##
runlist.runname = runname.(eachrow(runlist))
leftjoin!(runlist, DataFrame(translations), on=:runname)
##
PIVlist = loadmeta(m -> !ismissing(m.TSI_ID) && !ismissing(m.registration_path))
Vortex.add_registrations!(PIVlist)
PIVlist.runname = runname.(eachrow(PIVlist))
leftjoin!(PIVlist, DataFrame(translations), on=:runname)
transform!(PIVlist, :tx => ByRow(x -> ismissing(x) ? 0 : x) => :tx)
transform!(PIVlist, :ty => ByRow(x -> ismissing(x) ? 0 : x) => :ty)

foreach(eachrow(PIVlist)) do runmeta
    data, path = produce_or_load(runmeta, datadir("PIV", "registered");
                     filename=runname, tag=false) do runmeta
        cine_PIV = if ismissing(runmeta.cine_ID)
            zeros(Float32, 768, 576)
        else
            cine = Vortex.phantom_bgsub(runmeta)
            i_PIV, t_err = if !ismissing(runmeta.timings_path)
                Vortex.phantom_timing(runmeta, Vortex.read_timings(runmeta), size(cine, 3))
            else
                (1, 0)
            end
            t_err > 0 && @warn "PIV time is not simultaneous with cine frame"
            cine[:, :, i_PIV]
        end
        uv = Vortex.transform_prana(runmeta, axes(cine_PIV); tx=-runmeta.tx, ty=runmeta.ty)
        phantomscale = runmeta.phantomscale
        xx = range(0, step=phantomscale.mm_per_px, length=size(cine_PIV, 2))
        yy = range(stop = phantomscale.origin_mm, step=phantomscale.mm_per_px, length=size(cine_PIV, 1))
        u, v = (reverse.((uv.u, uv.v), dims=2)...,)
        data = Dict{String, Any}("x" => xx, "y" => yy, "u" => u, "v" => v)
    end
    data["runname"] = runname(runmeta)
    data
end

##
f = Figure(resolution=(1100,1300))

ax = Axis(f[1, 1], aspect=DataAspect())
cine_adjusted = rotr90(imadjust(Float32.(cine_PIV), qmax=0.9995))
phantomscale = runmeta.phantomscale
xx = range(0, step=phantomscale.mm_per_px, length=size(cine_adjusted, 1))
yy = range(stop = phantomscale.origin_mm, step=phantomscale.mm_per_px, length=size(cine_adjusted, 2))
heatmap!(ax, xx, yy, cine_adjusted, 
    colormap=[RGBA(1, 1, 1, 1), RGBA(0, 0, 0, 1)])

v_shocked = filter(<(-200), uv.v)
V₀ = isempty(v_shocked) ? 0 : median(v_shocked)
v′ = copy(uv.v)
v′[v′ .< V₀/2] .-= V₀
function u_xy(xx, yy, u, v)
    u_itp = linear_interpolation((xx, yy), u, extrapolation_bc=NaN)
    v_itp = linear_interpolation((xx, yy), v, extrapolation_bc=NaN)
    (x, y) -> Point2(u_itp(x, y), v_itp(x, y))
end
u = u_xy(xx, yy, reverse.((uv.u, v′), dims=2)...)

streamplot!(ax, u, xx, yy, stepsize=0.5, gridsize=(128, 128, 128), arrow_size=5)
f
