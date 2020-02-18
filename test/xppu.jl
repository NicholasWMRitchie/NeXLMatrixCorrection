using Revise
using NeXLCore
using NeXLMatrixCorrection
using Random
using Test

@testset "XPPU" begin
    rgen = MersenneTwister(0xBADF00D)
    material, cxr, e0 = "K240", n"O K-L3", 20.0e3
    k240_mat = NeXLCore.material("K240", #
                Dict(n"O"=>0.340023, n"Mg"=>0.030154, n"Si"=>0.186986, n"Ti"=>0.059950, #
                n"Zn"=>0.040168, n"Zr"=>0.074030, n"Ba"=>0.268689),missing)
    k240_af = atomicfraction(k240_mat)

    m=NeXLMatrixCorrection.m(inner(cxr))

    k240 =  uvs(
        ( MassFractionLabel(material, elm)=>uv(k240_mat[elm],0.01*k240_mat[elm]) for elm in keys(k240_mat))...,
        ( AtomicWeightLabel(material, elm)=>uv(a(elm, k240_mat),0.001*a(elm,k240_mat)) for elm in keys(k240_mat))...,
        ( JzLabel(elm)=> Ju(elm) for elm in keys(k240_mat) )...,
        NeXLMatrixCorrection.E0Label(material)=>uv(e0,0.1e3),
        NeXLMatrixCorrection.mLabel(inner(cxr))=>uv(m,0.01*m)
    );

    allinp = AllInputs()
    mjz = StepMJZbarb(material, [ n"O", n"Mg", n"Si", n"Ti", n"Zn", n"Zr", n"Ba" ] ) |
        MaintainLabels([NeXLMatrixCorrection.mLabel], k240)

    mjz_res = mjz(k240);
    mjz_mcres = mcpropagate(mjz, k240, 10000, parallel=false, rng=rgen);

    @test isapprox(value(mjz_res[NeXLMatrixCorrection.BigMLabel(material)]),0.466,atol=0.001)
    @test isapprox(value(mjz_res[NeXLMatrixCorrection.ZbarbLabel(material)]),22.5,atol=0.1)
    @test isapprox(value(mjz_res[NeXLMatrixCorrection.E0keVLabel(material)]),20.0, atol=0.1)
    @test isapprox(value(mjz_res[NeXLMatrixCorrection.JLabel(material)]),0.217, atol=0.001)
    @test isapprox(σ(mjz_res[NeXLMatrixCorrection.BigMLabel(material)]),0.0023,atol=0.0001)
    @test isapprox(σ(mjz_res[NeXLMatrixCorrection.ZbarbLabel(material)]),0.25,atol=0.03)
    @test isapprox(σ(mjz_res[NeXLMatrixCorrection.E0keVLabel(material)]),0.1, atol=0.01)
    @test isapprox(σ(mjz_res[NeXLMatrixCorrection.JLabel(material)]),0.00198, atol=0.0001)

    # println("Analytical Result")
    # print(mjz_res)
    # println("MC Result")
    # print(mjz_mcres)

    dpt = StepDPT(material, inner(cxr)) | allinp

    dpt_res = dpt(mjz_res);
    dpt_mcres = mcpropagate(dpt ∘ mjz, k240, 1000, parallel=false, rng=rgen);

    @test isapprox(value(dpt_res[NeXLMatrixCorrection.DkLabel(material,inner(cxr),1)]),6.6e-6,atol=0.01e-6)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.DkLabel(material,inner(cxr),2)]),1.49e-05,atol=0.01e-5)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.DkLabel(material,inner(cxr),3)]),1.02e-5,atol=0.01e-5)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.PkLabel(material,inner(cxr),1)]),0.780,atol=0.001)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.PkLabel(material,inner(cxr),2)]),0.100,atol=0.001)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.PkLabel(material,inner(cxr),3)]),-0.446,atol=0.001)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.TkLabel(material,inner(cxr),1)]),0.911,atol=0.001)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.TkLabel(material,inner(cxr),2)]),0.231,atol=0.001)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.TkLabel(material,inner(cxr),3)]),-0.315,atol=0.001)

    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.DkLabel(material,inner(cxr),1)]),0.0,atol=0.001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.DkLabel(material,inner(cxr),2)]),4.3e-09,atol=0.1e-9)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.DkLabel(material,inner(cxr),3)]),9.3e-08,atol=0.1e-8)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.PkLabel(material,inner(cxr),1)]),0.0,atol=0.001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.PkLabel(material,inner(cxr),2)]),0.0,atol=0.001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.PkLabel(material,inner(cxr),3)]),0.00049,atol=0.00001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.TkLabel(material,inner(cxr),1)]),0.009,atol=0.001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.TkLabel(material,inner(cxr),2)]),0.009,atol=0.001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.TkLabel(material,inner(cxr),3)]),0.0088,atol=0.0001)


    @test isapprox(value(dpt_res[NeXLMatrixCorrection.BigMLabel(material)]),0.466,atol=0.001)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.ZbarbLabel(material)]),22.5,atol=0.1)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.E0keVLabel(material)]),20.0, atol=0.1)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.JLabel(material)]),0.217, atol=0.001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.BigMLabel(material)]),0.0023,atol=0.0001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.ZbarbLabel(material)]),0.25,atol=0.03)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.E0keVLabel(material)]),0.1, atol=0.01)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.JLabel(material)]),0.00198, atol=0.0001)

    # println("Analytical Result")
    # print(dpt_res)
    # println("MC Result")
    # print(dpt_mcres)

    retain = MaintainLabels([NeXLMatrixCorrection.E0keVLabel, NeXLMatrixCorrection.ZbarbLabel ], dpt_res)

    qla = StepQlaOoS(material, inner(cxr)) | retain
    qla_res = qla(dpt_res);
    qla_model = qla ∘ dpt ∘ mjz
    qla_mcres = mcpropagate(qla_model, k240, 10000, parallel=false, rng=rgen)

    #println("Analytical Result")
    #print(qla_res)
    #println("MC Result")
    #print(qla_mcres)

    @test isapprox(value(qla_res[NeXLMatrixCorrection.QlaLabel(material,inner(cxr))]), 0.548, atol=0.001)
    @test isapprox(value(qla_res[NeXLMatrixCorrection.OoSLabel(material, inner(cxr))]), 0.00108, atol=0.00001)
    @test isapprox(value(qla_res[NeXLMatrixCorrection.ηLabel(material)]),0.253, atol=0.001)
    @test isapprox(value(qla_res[NeXLMatrixCorrection.JU0Label(material)]), 99.8, atol=0.1)
    @test isapprox(value(qla_res[NeXLMatrixCorrection.WbarLabel(material)]), 0.665, atol=0.001)
    @test isapprox(value(qla_res[NeXLMatrixCorrection.qLabel(material)]),0.987, atol=0.001)

    @test isapprox(σ(qla_res[NeXLMatrixCorrection.QlaLabel(material,inner(cxr))]), 0.0173, atol=0.001)
    @test isapprox(σ(qla_res[NeXLMatrixCorrection.OoSLabel(material, inner(cxr))]), 2.85e-5, atol=0.1e-5)
    @test isapprox(σ(qla_res[NeXLMatrixCorrection.ηLabel(material)]),0.00218, atol=0.0001)
    @test isapprox(σ(qla_res[NeXLMatrixCorrection.JU0Label(material)]), 0.682, atol=0.001)
    @test isapprox(σ(qla_res[NeXLMatrixCorrection.WbarLabel(material)]), 0.000667, atol=0.00001)
    @test isapprox(σ(qla_res[NeXLMatrixCorrection.qLabel(material)]),0.00595, atol=0.0001)

    maintain = MaintainLabels([NeXLMatrixCorrection.E0keVLabel, NeXLMatrixCorrection.ZbarbLabel, NeXLMatrixCorrection.OoSLabel, NeXLMatrixCorrection.QlaLabel ], qla_res)

    rp = NeXLMatrixCorrection.StepRPhi0(material, inner(cxr)) | maintain
    rp_res = rp(qla_res)
    rp_model = rp ∘ qla_model
    rp_mcres = mcpropagate(rp_model, k240, 1000, parallel=false, rng=rgen)

    #println("Analytical Result")
    #print(rp_res)
    #println("MC Result")
    #print(rp_mcres)

    xpp = NeXLMatrixCorrection.XPP(k240_mat,inner(cxr),e0)

    @test isapprox(value(rp_res[NeXLMatrixCorrection.RLabel(material,inner(cxr))]), 0.852, atol=0.001)
    @test isapprox(value(rp_res[NeXLMatrixCorrection.ϕ0Label(material,inner(cxr))]), xpp.ϕ0, atol=0.01)

    @test isapprox(σ(rp_res[NeXLMatrixCorrection.RLabel(material,inner(cxr))]), 0.00146, atol=0.001)
    @test isapprox(σ(rp_res[NeXLMatrixCorrection.ϕ0Label(material, inner(cxr))]), 0.00643, atol=0.0001)

    maintain = MaintainLabels([NeXLMatrixCorrection.ZbarbLabel, NeXLMatrixCorrection.ϕ0Label, NeXLMatrixCorrection.E0keVLabel], rp_res)

    frbar = NeXLMatrixCorrection.StepFRBar(material, inner(cxr)) | maintain
    frbar_res = frbar(rp_res)
    frbar_model = frbar ∘ rp_model
    frbar_mcres = mcpropagate(frbar_model, k240, 1000, parallel=false, rng=rgen)

    #println("Analytical Result")
    #print(frbar_res)
    #println("MC Result")
    #print(frbar_mcres)

    @test isapprox(value(frbar_res[NeXLMatrixCorrection.FLabel(material,inner(cxr))]), xpp.F, atol=0.01*xpp.F)
    @test isapprox(value(frbar_res[NeXLMatrixCorrection.RbarLabel(material,inner(cxr))]), 3.31e-4, atol=0.01e-4)
    @test isapprox(σ(frbar_res[NeXLMatrixCorrection.FLabel(material,inner(cxr))]), 2.56e-5, atol=0.1e-5)
    @test isapprox(σ(frbar_res[NeXLMatrixCorrection.RbarLabel(material, inner(cxr))]), 5.36e-6, atol=0.1e-6)

    maintain = MaintainLabels([NeXLMatrixCorrection.ϕ0Label, NeXLMatrixCorrection.RbarLabel, NeXLMatrixCorrection.FLabel], frbar_res)

    pb = NeXLMatrixCorrection.StepPb(material, inner(cxr)) | maintain
    pb_res = pb(frbar_res)
    pb_model = pb ∘ frbar_model
    pb_mcres = mcpropagate(pb_model, k240, 1000, parallel=false, rng=rgen)

    @test isapprox(value(pb_res[NeXLMatrixCorrection.PLabel(material,inner(cxr))]), 1.42e4, atol=0.01e4)
    @test isapprox(value(pb_res[NeXLMatrixCorrection.bLabel(material,inner(cxr))]), 7.78e3, atol=0.01e3)
    @test isapprox(σ(pb_res[NeXLMatrixCorrection.PLabel(material,inner(cxr))]), 2.85e2, atol=0.4e2)
    @test isapprox(σ(pb_res[NeXLMatrixCorrection.bLabel(material,inner(cxr))]), 1.25e2, atol=0.2e2)

    aϵ = NeXLMatrixCorrection.Stepaϵ(material, inner(cxr))
    aϵ_res = aϵ(pb_res)
    aϵ_model = aϵ ∘ pb_model
    aϵ_mcres = mcpropagate(aϵ_model, k240, 1000, parallel=false, rng=rgen)

    @test isapprox(value(aϵ_res[NeXLMatrixCorrection.aLabel(material,inner(cxr))]), 6.75e3, atol=0.01e3)
    @test isapprox(value(aϵ_res[NeXLMatrixCorrection.ϵLabel(material,inner(cxr))]), -0.133, atol=0.001)
    @test isapprox(σ(aϵ_res[NeXLMatrixCorrection.aLabel(material,inner(cxr))]), 2.85e2, atol=0.4e2)
    @test isapprox(σ(aϵ_res[NeXLMatrixCorrection.ϵLabel(material,inner(cxr))]), 1.25e2, atol=0.2e2)


    @test isapprox(value(A_res[NeXLMatrixCorrection.ALabel(material,inner(cxr))]), 393, atol=1.0)

end
