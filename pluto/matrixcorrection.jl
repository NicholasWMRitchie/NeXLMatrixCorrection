### A Pluto.jl notebook ###
# v0.12.7

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

# ╔═╡ 5040ab20-2050-11eb-2fc1-2decae24bed9
begin
	using NeXLMatrixCorrection
	using DataFrames, PlutoUI, CSV
	md"""
	## Matrix Correction Calculator
	"""
end

# ╔═╡ 65997380-2050-11eb-3e3d-b775f365f826
md"""
#### Sample and Standard Material

Enter the composition of the sample and standard.  Compositions may be entered as chemical formulae as in the default case.  Alternatively, compositions may be entered as mass-fractions like "0.54753\*Na+0.45247\*F".

Sample:   $(@bind sampstr TextField(default="NaAlSi3O8"))

Standard: $(@bind stdstr TextField(default="NaF"))
"""

# ╔═╡ e4353a00-2053-11eb-0380-03ab59b5f6b3
md"""
#### Instrument Parameters

Beam energy (keV) $(@bind e0 NumberField(1.0:0.1:100.0, default=20.0))

Take-off angle (°) $(@bind toa NumberField(0.0:0.1:90.0, default=40.0))
"""

# ╔═╡ 46ccada0-2055-11eb-258e-251301e807c6
begin
	sample = parse(Material, sampstr)
	standard = parse(Material, stdstr)
	md"""
	#### Sample
	$(asa(DataFrame, sample))
	#### Standard
	$(asa(DataFrame, standard))
	"""
end

# ╔═╡ c32171b0-2055-11eb-299b-fd7ea9d83b83
begin
	els=intersect(keys(standard), keys(sample))
	dfs=DataFrame[]
	for elm in els
		for trs in (ktransitions, ltransitions, mtransitions)
			cxrs = characteristic(elm, trs, 0.01, 1.0e3*e0)
			if length(cxrs)>0
				zafs = zafcorrection(
					XPP, ReedFluorescence, NullCoating,
					sample, standard, cxrs, 1.0e3*e0)
				push!(dfs,asa(DataFrame, zafs..., deg2rad(toa), deg2rad(toa)))

			end
		end
	end
	edf=vcat(dfs...)
	ec=CSV.write(joinpath(tempdir(),"EDS.csv"), edf)
	md"""
	#### EDS Mode Results
	$(edf[:,5:end])
	Written to: $ec
	"""
end

# ╔═╡ f11e5580-2058-11eb-084c-1dea7dcd01e0
begin
	wdf=asa(DataFrame, 
		sample, Dict(el=>standard for el in els),
		e0*1.0e3, deg2rad(toa), deg2rad(toa), 
		mctype=XPP, fctype=ReedFluorescence, cctype=Coating)
	cw=CSV.write(joinpath(tempdir(),"WDS.csv"), wdf)
	md"""
	#### WDS Mode Results
	
	$(wdf[:,7:end])
	Written to: $cw
	"""
end

# ╔═╡ Cell order:
# ╟─5040ab20-2050-11eb-2fc1-2decae24bed9
# ╟─65997380-2050-11eb-3e3d-b775f365f826
# ╟─e4353a00-2053-11eb-0380-03ab59b5f6b3
# ╟─46ccada0-2055-11eb-258e-251301e807c6
# ╟─c32171b0-2055-11eb-299b-fd7ea9d83b83
# ╟─f11e5580-2058-11eb-084c-1dea7dcd01e0
