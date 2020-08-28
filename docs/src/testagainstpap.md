### Testing Matrix Correction Algorithms
Testing matrix correction algorithms XPP and CitZAF against the Pouchou and Pichoir k-ratio database.

````julia
using CSV
using DataFrames
using NeXLMatrixCorrection

function mapline(elm, index)
  trss = Dict(0=>kalpha, 2=>kbeta, 12=>lalpha, 31=>lbeta, 72=>malpha, 69=>mbeta)
  trs=trss[index]
  return [ brightest(characteristic(elm, trs)) ]
end

pap = CSV.read("papkratios.csv", DataFrame, header=3, skipto=4)
xppres, czres = Union{Float64,Missing}[], Union{Float64,Missing}[]
for r in eachrow(pap)
  try
    a, b = elements[r.A], elements[r.B]
    e0, θ  = 1.0e3*r.E0, deg2rad(r.TOA)
    std, unk = pure(a), material("Unknown",Dict(a=>r.WgtFracA, b=>1.0-r.WgtFracA))
    kk, lines = r.kA, mapline(a, r.Line)
    algs = zafcorrection(XPP, ReedFluorescence, NullCoating, unk, std, lines, e0)
    push!(xppres, k(algs..., θ, θ)/kk)
    algs = zafcorrection(CitZAF, ReedFluorescence, NullCoating, unk, std, lines, e0)
    push!(czres, k(algs..., θ, θ)/kk)
  catch
    push!(xppres, missing)
    push!(czres, missing)
  end
end
insertcols!(pap, ncol(pap)+1, :XPP=>xppres)
insertcols!(pap, ncol(pap)+1, :CitZAF=>czres)
display(pap)
````


````
826×9 DataFrame. Omitted printing of 2 columns
│ Row │ A     │ Line  │ B     │ E0      │ WgtFracA │ kA      │ TOA     │
│     │ Int64 │ Int64 │ Int64 │ Float64 │ Float64  │ Float64 │ Float64 │
├─────┼───────┼───────┼───────┼─────────┼──────────┼─────────┼─────────┤
│ 1   │ 13    │ 0     │ 26    │ 20.0    │ 0.241    │ 0.124   │ 52.5    │
│ 2   │ 13    │ 0     │ 26    │ 25.0    │ 0.241    │ 0.098   │ 52.5    │
│ 3   │ 13    │ 0     │ 26    │ 30.0    │ 0.241    │ 0.083   │ 52.5    │
│ 4   │ 26    │ 0     │ 13    │ 20.0    │ 0.759    │ 0.736   │ 52.5    │
│ 5   │ 26    │ 0     │ 13    │ 25.0    │ 0.759    │ 0.742   │ 52.5    │
│ 6   │ 26    │ 0     │ 13    │ 30.0    │ 0.759    │ 0.748   │ 52.5    │
│ 7   │ 26    │ 0     │ 16    │ 10.0    │ 0.466    │ 0.406   │ 75.0    │
⋮
│ 819 │ 42    │ 12    │ 7     │ 6.0     │ 0.942    │ 0.9023  │ 40.0    │
│ 820 │ 42    │ 12    │ 7     │ 8.0     │ 0.942    │ 0.9078  │ 40.0    │
│ 821 │ 42    │ 12    │ 7     │ 10.0    │ 0.942    │ 0.9122  │ 40.0    │
│ 822 │ 42    │ 12    │ 7     │ 12.0    │ 0.942    │ 0.9163  │ 40.0    │
│ 823 │ 42    │ 12    │ 7     │ 15.0    │ 0.942    │ 0.922   │ 40.0    │
│ 824 │ 42    │ 12    │ 7     │ 20.0    │ 0.942    │ 0.9268  │ 40.0    │
│ 825 │ 42    │ 12    │ 7     │ 25.0    │ 0.942    │ 0.936   │ 40.0    │
│ 826 │ 42    │ 12    │ 7     │ 30.0    │ 0.942    │ 0.9399  │ 40.0    │
````





##### XPP
Let's visualize this.
````julia
using Gadfly
plot(pap, x=:XPP, y=:XPP, Stat.histogram(bincount=50), Geom.bar, Guide.title("XPP"))
````


![](figures/testagainstpap_2_1.svg)



##### CitZAF
````julia
plot(pap, x=:CitZAF, y=:CitZAF, Stat.histogram(bincount=50), Geom.bar, Guide.title("CitZAF"))
````


![](figures/testagainstpap_3_1.svg)



##### Summary Statistics
````julia
describe(pap[:,end-1:end], :mean, :std, :min, :q25, :median, :q75, :max)
````


````
2×8 DataFrame. Omitted printing of 2 columns
│ Row │ variable │ mean     │ std       │ min      │ q25      │ median   │
│     │ Symbol   │ Float64  │ Float64   │ Float64  │ Float64  │ Float64  │
├─────┼──────────┼──────────┼───────────┼──────────┼──────────┼──────────┤
│ 1   │ XPP      │ 1.00355  │ 0.0264633 │ 0.897411 │ 0.989226 │ 0.999223 │
│ 2   │ CitZAF   │ 0.990101 │ 0.0477465 │ 0.848773 │ 0.96709  │ 0.995988 │
````
