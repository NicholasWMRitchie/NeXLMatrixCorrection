### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 6538daa0-243a-11eb-2726-c9068913ab03
begin
	using PlutoUI, NeXLMatrixCorrection, Gadfly
	md"""
	# Plot the ϕ(ρz) Curve
	"""
end

# ╔═╡ a3e07a00-243b-11eb-01e9-27ea5be6f71e
md"""
## Measurement Parameters

Material: $(@bind matstr TextField(default="Al"))

Beam energy (keV): $(@bind e0keV NumberField(0.0:0.1:100.0,default=29.0))

Characteristic X-ray: $(@bind cxrstr TextField(default="Mg K-L3"))

Take-Off Angle (°): $(@bind toa NumberField(1.0:0.1:89.0, default=40.0))
"""

# ╔═╡ 8f049f40-243a-11eb-2424-ef70b99f4423
plot([ XPP, CitZAF, Riveros1993, XPhi ], parse(Material,matstr), parse(CharXRay, cxrstr), 1.0e3*e0keV, deg2rad(toa))

# ╔═╡ c0ad8e42-2826-11eb-1b81-9f9a11b6c51d
begin
	cxr = parse(CharXRay,cxrstr)
	sh = inner(cxr)
md"""
  * E($cxr) = $(0.001*energy(cxr)) keV
  * Ec($sh) = $(0.001*energy(sh)) keV
  * U($sh @ $e0keV keV) = $(1.0e3*e0keV/energy(sh))
"""
end


# ╔═╡ Cell order:
# ╟─6538daa0-243a-11eb-2726-c9068913ab03
# ╟─a3e07a00-243b-11eb-01e9-27ea5be6f71e
# ╟─8f049f40-243a-11eb-2424-ef70b99f4423
# ╟─c0ad8e42-2826-11eb-1b81-9f9a11b6c51d
