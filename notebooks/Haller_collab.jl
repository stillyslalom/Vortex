### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ cd2395ce-d435-11ed-26a0-75f24307f5cc
begin
	using Pkg; Pkg.activate("../.")
	using Vortex
	using PlutoUI
	using ImageCore
	using ImageShow
end

# ╔═╡ 36e893ea-3674-4331-84ba-8b2adbdf9867
runmeta = loadmeta() do m
	# Load initial condition runs with PLIF data
	!ismissing(m.cine_ID) && !ismissing(m.timings_path) && ismissing(m.ptrace_path)
end;

# ╔═╡ 1678d853-ab70-4aeb-9de2-89186055f0e4
@bind gas Select(["N2", "Ar", "CF4", "SF6"])

# ╔═╡ ce7a1da1-c783-47d3-8216-473a4d973aba
matchruns = filter(m -> occursin(gas, m.MST_gas), runmeta)

# ╔═╡ 21c55f66-f492-4946-a177-053961964990
@bind run Select([string(r.Date)*'_'*r.ID for r in eachrow(matchruns)])

# ╔═╡ 325d8965-bf7e-4c4f-86a4-edfa61ccd682
cine = phantom_bgsub(select_run(runmeta, run));

# ╔═╡ ca78e0a6-693a-40cf-bcae-ba979a2b95e6
@bind i Clock(0.05, false, false, size(cine, 3))

# ╔═╡ c8bcba1b-df29-4319-a0e6-e73530418b7c
"t = $(round(i*0.05, sigdigits=5)) ms"

# ╔═╡ 3084480e-e20b-45e5-a650-8b9497e6b28b
reinterpret(Gray{Float64}, imadjust(@view(cine[:,:,i]), qmax=0.9995))

# ╔═╡ Cell order:
# ╟─cd2395ce-d435-11ed-26a0-75f24307f5cc
# ╠═36e893ea-3674-4331-84ba-8b2adbdf9867
# ╠═ce7a1da1-c783-47d3-8216-473a4d973aba
# ╠═1678d853-ab70-4aeb-9de2-89186055f0e4
# ╠═21c55f66-f492-4946-a177-053961964990
# ╠═325d8965-bf7e-4c4f-86a4-edfa61ccd682
# ╟─ca78e0a6-693a-40cf-bcae-ba979a2b95e6
# ╟─c8bcba1b-df29-4319-a0e6-e73530418b7c
# ╠═3084480e-e20b-45e5-a650-8b9497e6b28b
