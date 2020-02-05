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
    MaintainLabels([NeXLMatrixCorrection.mLabel, NeXLMatrixCorrection.E0Label], k240)

mjz_res = mjz(k240);
mjz_mcres = mcpropagate(mjz, k240, 10000, parallel=false, rng=rgen);

println("Analytical Result")
print(mjz_res)
println("MC Result")
print(mjz_mcres)

dpt = StepDPT(material, inner(cxr)) | allinp

dpt_res = dpt(mjz_res);
dpt_mcres = mcpropagate(dpt, mjz_res, 1000, parallel=false, rng=rgen);

println("Analytical Result")
print(dpt_res)
println("MC Result")
print(dpt_mcres)

qla = StepQlaOoS(material, inner(cxr))
qla_res = qla(dpt_res);
qla_mcres = mcpropagate(qla, dpt_res, 1000, parallel=false, rng=rgen)

println("Analytical Result")
print(qla_res)
println("MC Result")
print(qla_mcres)




#end
