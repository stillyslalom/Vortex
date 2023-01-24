### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 9b52f6d2-8b7a-11ed-37a5-3d6cf4bd72bd
begin
	using Pkg
	Pkg.activate("../.")
	using PressureTraceXT: PressureTrace, xt, shockspeed
	using DataFramesMeta
	using Interpolations
	using CairoMakie
	using Measurements
	using Unitful
	using Statistics
	using PyThermo
	using PyThermo.ShockTube
end

# ╔═╡ 6665742a-e1d9-41cc-9495-1533b5826447
begin
	# ptrace_path = raw"S:\users\shocktube\runs\2022\11\17\run2\ptrace.lvm"
	ptrace_path = raw"S:\users\shocktube\runs\2023\01\11\run4\ptrace.lvm"
	# ptrace_path = raw"S:\users\shocktube\runs\2022\9\12\rupture_18ga_steel.lvm"
	# ptrace_path = raw"S:\users\shocktube\runs\2023\01\17\18ga_rupture.lvm"

	channelmap_path =
	raw"S:\users\shocktube\hardware_laser_camera_windows_optics_etc\PT_locations\channelmaps\2022-12-08.csv"
	# channelmap_path = 
	# 	raw"S:\users\shocktube\hardware_laser_camera_windows_optics_etc\PT_locations\channelmaps\2021-10-08.csv"
end

# ╔═╡ f7f9a163-cefd-4698-868a-92c1c77c8527
begin
	ptrace = PressureTrace(ptrace_path, channelmap_path)
	select!(ptrace.data, Not(:PT9))
end

# ╔═╡ a2dfb7ae-7e91-4800-afb7-1ab3dd9d4613
xtdata = xt(ptrace,
raw"S:\users\shocktube\hardware_laser_camera_windows_optics_etc\PT_locations\2022-09-12.csv", :PT3 => 20)

# ╔═╡ 338d0482-1080-4d27-8144-853ef09bbc5d
plot(Measurements.value.(xtdata.x), xtdata.t)

# ╔═╡ 2f7cb87a-0981-47bb-a19f-c3f8e698c795
itp = LinearInterpolation(xtdata.x, xtdata.t)

# ╔═╡ 38dbf891-055a-47bf-8688-3d2a3e459e9c
begin
	1e6*(itp(xtdata.x[end] - 0.295) - xtdata.t[2])
end

# ╔═╡ b11de9c0-6812-4b62-9ad4-8309cd1dfbbd
29.5u"cm" |> u"inch"

# ╔═╡ ecc85df2-c68b-41a9-93ad-90c0b3421798
driver = Species("N2", P=14.4u"psi", T=21u"°C")

# ╔═╡ d10a7bac-ce3f-4906-a259-54a328f38555
driven = Species("N2", P=14.4u"psi", T=21u"°C")

# ╔═╡ 1bcae8e1-4ae8-45ab-99fa-f53e832db399
Wₛ = Measurements.value(shockspeed(xtdata)[3:end] |> mean)

# ╔═╡ edd859f8-ddd4-40a5-bec2-59a8f2067cbe
sc = shockcalc(driver, driven, ustrip(Wₛ / soundspeed(driven)))

# ╔═╡ 3806bf98-7d74-4a4f-8585-279a362f9fa6
md"""
Account for driver-to-driven constriction - WiSTL calculator uses a factor of 1.35, but 1.2 seems to fit data better
"""

# ╔═╡ 88476cb0-a62e-45af-acff-d3c790f6227d
p_driver = (pressure(sc.driver)/1.2 - 14.5u"psi") |> u"psi"

# ╔═╡ Cell order:
# ╠═9b52f6d2-8b7a-11ed-37a5-3d6cf4bd72bd
# ╠═f7f9a163-cefd-4698-868a-92c1c77c8527
# ╠═a2dfb7ae-7e91-4800-afb7-1ab3dd9d4613
# ╠═338d0482-1080-4d27-8144-853ef09bbc5d
# ╠═6665742a-e1d9-41cc-9495-1533b5826447
# ╠═38dbf891-055a-47bf-8688-3d2a3e459e9c
# ╠═2f7cb87a-0981-47bb-a19f-c3f8e698c795
# ╠═b11de9c0-6812-4b62-9ad4-8309cd1dfbbd
# ╠═ecc85df2-c68b-41a9-93ad-90c0b3421798
# ╠═d10a7bac-ce3f-4906-a259-54a328f38555
# ╠═1bcae8e1-4ae8-45ab-99fa-f53e832db399
# ╠═edd859f8-ddd4-40a5-bec2-59a8f2067cbe
# ╟─3806bf98-7d74-4a4f-8585-279a362f9fa6
# ╠═88476cb0-a62e-45af-acff-d3c790f6227d
