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

println("Analytical Result")
print(mjz_res)
println("MC Result")
print(mjz_mcres)

dpt = StepDPT(material, inner(cxr)) | allinp

dpt_res = dpt(mjz_res);
dpt_mcres = mcpropagate(dpt ∘ mjz, k240, 1000, parallel=false, rng=rgen);

println("Analytical Result")
print(dpt_res)
println("MC Result")
print(dpt_mcres)

qla = StepQlaOoS(material, inner(cxr))
qla_res = qla(dpt_res);
qla_model = qla ∘ dpt ∘ mjz
qla_mcres = mcpropagate(qla_model, k240, 1000, parallel=false, rng=rgen)

println("Analytical Result")
print(qla_res)
println("MC Result")
print(qla_mcres)

struct StepRPhi0 <: MeasurementModel
    material::String
    shell::AtomicSubShell
end

struct RLabel <: Label
    material::String
end

struct ϕ0Label <: Label
    material::String
    shell::AtomicSubShell
end

function NeXLUncertainties.compute(rp::StepRPhi0, inputs::LabeledValues, withJac::Bool)::MMResult
    # inputs
    E0l, ql, JU0l = E0keVLabel(rp.material), qLabel(rp.Material), JU0Label(rp.material)
    ηl, Wbarl = ηLabel(rp.Material), WbarLabel(rp.Material)
    e0, q, JU0 = inputs[E0l], inputs[ql], inputs[JU0l]
    η, Wbar = inputs[ηl], inputs[Wbarl]
    Ea = 0.001 * energy(rp.shell)
    U0 = e0/Ea
    # outputs
    GU0 = (U0-1.0-(1.0-1.0/U0^(1.0+q))/(1+q))/((2.0+q)*JU0)
    R = 1.0 - η*Wbar*(1.0-GU0)
    ϕ0 = 1.0+3.3*(1.0-1.0/U0^(2.0-2.3*η))*η^1.2

    Rl, ϕ0l = RLabel(rp.material), ϕ0Label(rp.material, rp.shell)
    vals = LabeledValues( [ Rl, ϕ0l ], [ R, ϕ0 ])
    jac = withJac ? zeros(Float64, length(vals), length(inputs)) : nothing
    if withJac
        δGUδE0 = (1.0/(JU0*Ea))*(((1.0-U0^(-2.0-q))/(2.0+q))-GU0)
        δGUδq = ((U0^(1.0+q)-1.0)-(1.0-q)*log(U0))/(JU0*(1.0+q)^2*(2.0+q)*U0^(1.0+q))
        jac[1, indexin(ηl, inputs)] = Wbar*(1.0-GU0) # δRδη
        jac[1, indexin(E0l, inputs)] = η*Wbar*δGUδE0 # δRδE0
        jac[1, indexin(ql, inputs)] = η*Wbar*δGUδq   # δRδq
        jac[1, indexin(Wbarl, inputs)] = η*(1.0-GU0) # δRδWbar
        jac[2, indexin(ηl, inputs)]= η^0.2*(3.96*(1.0-U0^(2.3η-2.0))-7.59*η*U0^(2.3η-2.0)*log(U0)) # δϕ0δη
        jac[2, indexin(E0l, inputs)] = 7.59*η^1.2*(0.869565-η)*U0^(2.3η-3.0) # δϕ0δE0
    end
    return (vals, jac)
end


#end
