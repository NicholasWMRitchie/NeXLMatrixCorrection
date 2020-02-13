using Revise
using NeXLCore
using NeXLMatrixCorrection
using Random
using Test

#@testset "K240" begin # Test xppu against xpp using K240

rgen = MersenneTwister(0xBADF00D)
material, cxr = "K240", n"O K-L3"
m=NeXLMatrixCorrection.m(inner(cxr))

k240 =  uvs(
    MassFractionLabel(material, n"O")=>uv(0.340023,0.0034), #
    MassFractionLabel(material, n"Mg")=>uv(0.030154,0.003), #
    MassFractionLabel(material, n"Si")=>uv(0.186986, 0.0019), #
    MassFractionLabel(material, n"Ti")=>uv(0.059950,0.001), #
    MassFractionLabel(material, n"Zn")=>uv(0.040168,0.001), #
    MassFractionLabel(material, n"Zr")=>uv(0.074030,0.001), #
    MassFractionLabel(material, n"Ba")=>uv(0.268689,0.0026), #

    AtomicWeightLabel(material, n"O")=>uv(0.340023,0.0034), #
    AtomicWeightLabel(material, n"Mg")=>uv(0.030154,0.003), #
    AtomicWeightLabel(material, n"Si")=>uv(0.186986, 0.0019), #
    AtomicWeightLabel(material, n"Ti")=>uv(0.059950,0.001), #
    AtomicWeightLabel(material, n"Zn")=>uv(0.040168,0.001), #
    AtomicWeightLabel(material, n"Zr")=>uv(0.074030,0.001), #
    AtomicWeightLabel(material, n"Ba")=>uv(0.268689,0.0026), #

    JzLabel(n"O") => Ju(n"O"), #
    JzLabel(n"Mg") => Ju(n"Mg"), #
    JzLabel(n"Si") => Ju(n"Si"), #
    JzLabel(n"Ti") => Ju(n"Ti"), #
    JzLabel(n"Zn") => Ju(n"Zn"), #
    JzLabel(n"Zr") => Ju(n"Zr"), #
    JzLabel(n"Ba") => Ju(n"Ba"), #

    NeXLMatrixCorrection.E0Label(material)=>uv(20.0e3,0.1e3),
    NeXLMatrixCorrection.mLabel(inner(cxr))=>uv(m,0.01*m)
);

allinp = AllInputs()
mjz = StepMJZbarb(material, [ n"O", n"Mg", n"Si", n"Ti", n"Zn", n"Zr", n"Ba" ] ) |
    MaintainLabels([NeXLMatrixCorrection.mLabel], k240)

mjz_res = mjz(k240);
mjz_mcres = mcpropagate(mjz, k240, 10000, parallel=false, rng=rgen);

@test isapprox(value(mjz_res[NeXLMatrixCorrection.BigMLabel(material)]),182.,atol=1.0)
@test isapprox(value(mjz_res[NeXLMatrixCorrection.ZbarbLabel(material)]),22.5,atol=0.1)
@test isapprox(value(mjz_res[NeXLMatrixCorrection.E0keVLabel(material)]),20.0, atol=0.1)
@test isapprox(value(mjz_res[NeXLMatrixCorrection.JLabel(material)]),0.339, atol=0.001)
@test isapprox(σ(mjz_res[NeXLMatrixCorrection.BigMLabel(material)]),2.3,atol=0.1)
@test isapprox(σ(mjz_res[NeXLMatrixCorrection.ZbarbLabel(material)]),0.25,atol=0.03)
@test isapprox(σ(mjz_res[NeXLMatrixCorrection.E0keVLabel(material)]),0.1, atol=0.01)
@test isapprox(σ(mjz_res[NeXLMatrixCorrection.JLabel(material)]),0.003, atol=0.001)

# println("Analytical Result")
# print(mjz_res)
# println("MC Result")
# print(mjz_mcres)

dpt = StepDPT(material, inner(cxr)) | allinp

dpt_res = dpt(mjz_res);
dpt_mcres = mcpropagate(dpt ∘ mjz, k240, 1000, parallel=false, rng=rgen);

@test isapprox(value(dpt_res[NeXLMatrixCorrection.DLabel(material,inner(cxr),1)]),6.6e-6,atol=0.01e-6)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.DLabel(material,inner(cxr),2)]),1.45e-05,atol=1.0e-7)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.DLabel(material,inner(cxr),3)]),6.5e-06,atol=0.1e-6)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.PLabel(material,inner(cxr),1)]),0.78,atol=0.01)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.PLabel(material,inner(cxr),2)]),0.1,atol=0.01)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.PLabel(material,inner(cxr),3)]),-0.415,atol=0.001)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.TLabel(material,inner(cxr),1)]),0.911,atol=0.001)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.TLabel(material,inner(cxr),2)]),0.231,atol=0.001)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.TLabel(material,inner(cxr),3)]),-0.285,atol=0.001)

@test isapprox(σ(dpt_res[NeXLMatrixCorrection.DLabel(material,inner(cxr),1)]),0.0,atol=0.001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.DLabel(material,inner(cxr),2)]),1.0e-08,atol=1.0e-9)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.DLabel(material,inner(cxr),3)]),6.0e-08,atol=1e-8)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.PLabel(material,inner(cxr),1)]),0.0,atol=0.001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.PLabel(material,inner(cxr),2)]),0.0,atol=0.001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.PLabel(material,inner(cxr),3)]),0.0007,atol=0.0001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.TLabel(material,inner(cxr),1)]),0.009,atol=0.001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.TLabel(material,inner(cxr),2)]),0.009,atol=0.001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.TLabel(material,inner(cxr),3)]),0.0088,atol=0.0001)

@test isapprox(value(dpt_res[NeXLMatrixCorrection.BigMLabel(material)]),182.,atol=1.0)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.ZbarbLabel(material)]),22.5,atol=0.1)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.E0keVLabel(material)]),20.0, atol=0.1)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.JLabel(material)]),0.339, atol=0.001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.BigMLabel(material)]),2.3,atol=0.1)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.ZbarbLabel(material)]),0.25,atol=0.03)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.E0keVLabel(material)]),0.1, atol=0.01)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.JLabel(material)]),0.003, atol=0.001)

# println("Analytical Result")
# print(dpt_res)
# println("MC Result")
# print(dpt_mcres)

retain = MaintainLabels([NeXLMatrixCorrection.E0keVLabel], dpt_res)

qla = StepQlaOoS(material, inner(cxr)) | retain
qla_res = qla(dpt_res);
qla_model = qla ∘ dpt ∘ mjz
qla_mcres = mcpropagate(qla_model, k240, 10000, parallel=false, rng=rgen)

println("Analytical Result")
print(qla_res)
println("MC Result")
print(qla_mcres)

rp = NeXLMatrixCorrection.StepRPhi0(material, inner(cxr))
rp_res = rp(qla_res)
rp_model = rp ∘ qla ∘ dpt ∘ mjz
rp_mcres = mcpropagate(rp_model, k240, 1000, parallel=false, rng=rgen)

#end
