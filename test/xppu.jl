using NeXLCore
using NeXLMatrixCorrection
using Random
using Test

#@testset "XPPU" begin
rgen = MersenneTwister(0xBADF00D)
matname, cxr, e0 = "K240", n"O K-L3", 20.0e3
k240_mat = NeXLCore.material(matname, #
            Dict(n"O"=>0.340023, n"Mg"=>0.030154, n"Si"=>0.186986, n"Ti"=>0.059950, #
            n"Zn"=>0.040168, n"Zr"=>0.074030, n"Ba"=>0.268689),missing)
k240_af = atomicfraction(k240_mat)

m=NeXLMatrixCorrection.m(inner(cxr))

k240 =  uvs(
    ( MassFractionLabel(matname, elm)=>uv(k240_mat[elm],0.01*k240_mat[elm]) for elm in keys(k240_mat))...,
    ( AtomicWeightLabel(matname, elm)=>uv(a(elm, k240_mat),0.001*a(elm,k240_mat)) for elm in keys(k240_mat))...,
    ( JzLabel(elm)=> Ju(elm) for elm in keys(k240_mat) )...,
    NeXLMatrixCorrection.E0Label(matname)=>uv(e0,0.1e3),
    NeXLMatrixCorrection.mLabel(inner(cxr))=>uv(m,0.01*m)
);

allinp = AllInputs()
mjz = StepMJZbarb(matname, [ n"O", n"Mg", n"Si", n"Ti", n"Zn", n"Zr", n"Ba" ] ) |
    MaintainLabels([NeXLMatrixCorrection.mLabel], k240)

mjz_res = mjz(k240);
mjz_mcres = mcpropagate(mjz, k240, 10000, parallel=false, rng=rgen);

@test isapprox(value(mjz_res[NeXLMatrixCorrection.BigMLabel(matname)]),0.466,atol=0.001)
@test isapprox(value(mjz_res[NeXLMatrixCorrection.ZbarbLabel(matname)]),22.5,atol=0.1)
@test isapprox(value(mjz_res[NeXLMatrixCorrection.E0keVLabel(matname)]),20.0, atol=0.1)
@test isapprox(value(mjz_res[NeXLMatrixCorrection.JLabel(matname)]),0.217, atol=0.001)
@test isapprox(σ(mjz_res[NeXLMatrixCorrection.BigMLabel(matname)]),0.0023,atol=0.0001)
@test isapprox(σ(mjz_res[NeXLMatrixCorrection.ZbarbLabel(matname)]),0.25,atol=0.03)
@test isapprox(σ(mjz_res[NeXLMatrixCorrection.E0keVLabel(matname)]),0.1, atol=0.01)
@test isapprox(σ(mjz_res[NeXLMatrixCorrection.JLabel(matname)]),0.00198, atol=0.0001)

# println("Analytical Result")
# print(mjz_res)
# println("MC Result")
# print(mjz_mcres)

dpt = StepDPT(matname, inner(cxr)) | allinp

dpt_res = dpt(mjz_res);
dpt_mcres = mcpropagate(dpt ∘ mjz, k240, 1000, parallel=false, rng=rgen);

@test isapprox(value(dpt_res[NeXLMatrixCorrection.DkLabel(matname,inner(cxr),1)]),6.6e-6,atol=0.01e-6)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.DkLabel(matname,inner(cxr),2)]),1.49e-05,atol=0.01e-5)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.DkLabel(matname,inner(cxr),3)]),1.02e-5,atol=0.01e-5)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.PkLabel(matname,inner(cxr),1)]),0.780,atol=0.001)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.PkLabel(matname,inner(cxr),2)]),0.100,atol=0.001)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.PkLabel(matname,inner(cxr),3)]),-0.446,atol=0.001)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.TkLabel(matname,inner(cxr),1)]),0.911,atol=0.001)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.TkLabel(matname,inner(cxr),2)]),0.231,atol=0.001)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.TkLabel(matname,inner(cxr),3)]),-0.315,atol=0.001)

@test isapprox(σ(dpt_res[NeXLMatrixCorrection.DkLabel(matname,inner(cxr),1)]),0.0,atol=0.001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.DkLabel(matname,inner(cxr),2)]),4.3e-09,atol=0.1e-9)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.DkLabel(matname,inner(cxr),3)]),9.3e-08,atol=0.1e-8)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.PkLabel(matname,inner(cxr),1)]),0.0,atol=0.001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.PkLabel(matname,inner(cxr),2)]),0.0,atol=0.001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.PkLabel(matname,inner(cxr),3)]),0.00049,atol=0.00001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.TkLabel(matname,inner(cxr),1)]),0.009,atol=0.001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.TkLabel(matname,inner(cxr),2)]),0.009,atol=0.001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.TkLabel(matname,inner(cxr),3)]),0.0088,atol=0.0001)


@test isapprox(value(dpt_res[NeXLMatrixCorrection.BigMLabel(matname)]),0.466,atol=0.001)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.ZbarbLabel(matname)]),22.5,atol=0.1)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.E0keVLabel(matname)]),20.0, atol=0.1)
@test isapprox(value(dpt_res[NeXLMatrixCorrection.JLabel(matname)]),0.217, atol=0.001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.BigMLabel(matname)]),0.0023,atol=0.0001)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.ZbarbLabel(matname)]),0.25,atol=0.03)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.E0keVLabel(matname)]),0.1, atol=0.01)
@test isapprox(σ(dpt_res[NeXLMatrixCorrection.JLabel(matname)]),0.00198, atol=0.0001)

# println("Analytical Result")
# print(dpt_res)
# println("MC Result")
# print(dpt_mcres)

retain = MaintainLabels([NeXLMatrixCorrection.E0keVLabel, NeXLMatrixCorrection.ZbarbLabel ], dpt_res)

qla = StepQlaOoS(matname, inner(cxr)) | retain
qla_res = qla(dpt_res);
qla_model = qla ∘ dpt ∘ mjz
qla_mcres = mcpropagate(qla_model, k240, 10000, parallel=false, rng=rgen)

#println("Analytical Result")
#print(qla_res)
#println("MC Result")
#print(qla_mcres)

@test isapprox(value(qla_res[NeXLMatrixCorrection.QlaLabel(matname,inner(cxr))]), 0.548, atol=0.001)
@test isapprox(value(qla_res[NeXLMatrixCorrection.OoSLabel(matname, inner(cxr))]), 0.00108, atol=0.00001)
@test isapprox(value(qla_res[NeXLMatrixCorrection.ηLabel(matname)]),0.253, atol=0.001)
@test isapprox(value(qla_res[NeXLMatrixCorrection.JU0Label(matname)]), 99.8, atol=0.1)
@test isapprox(value(qla_res[NeXLMatrixCorrection.WbarLabel(matname)]), 0.665, atol=0.001)
@test isapprox(value(qla_res[NeXLMatrixCorrection.qLabel(matname)]),0.987, atol=0.001)

@test isapprox(σ(qla_res[NeXLMatrixCorrection.QlaLabel(matname,inner(cxr))]), 0.0173, atol=0.001)
@test isapprox(σ(qla_res[NeXLMatrixCorrection.OoSLabel(matname, inner(cxr))]), 2.85e-5, atol=0.1e-5)
@test isapprox(σ(qla_res[NeXLMatrixCorrection.ηLabel(matname)]),0.00218, atol=0.0001)
@test isapprox(σ(qla_res[NeXLMatrixCorrection.JU0Label(matname)]), 0.682, atol=0.001)
@test isapprox(σ(qla_res[NeXLMatrixCorrection.WbarLabel(matname)]), 0.000667, atol=0.00001)
@test isapprox(σ(qla_res[NeXLMatrixCorrection.qLabel(matname)]),0.00595, atol=0.0001)

maintain = MaintainLabels([NeXLMatrixCorrection.E0keVLabel, NeXLMatrixCorrection.ZbarbLabel, NeXLMatrixCorrection.OoSLabel, NeXLMatrixCorrection.QlaLabel ], qla_res)

rp = NeXLMatrixCorrection.StepRPhi0(matname, inner(cxr)) | maintain
rp_res = rp(qla_res)
rp_model = rp ∘ qla_model
rp_mcres = mcpropagate(rp_model, k240, 1000, parallel=false, rng=rgen)

#println("Analytical Result")
#print(rp_res)
#println("MC Result")
#print(rp_mcres)

xpp = NeXLMatrixCorrection.XPP(k240_mat,inner(cxr),e0)

@test isapprox(value(rp_res[NeXLMatrixCorrection.RLabel(matname,inner(cxr))]), 0.852, atol=0.001)
@test isapprox(value(rp_res[NeXLMatrixCorrection.ϕ0Label(matname,inner(cxr))]), xpp.ϕ0, atol=0.01)

@test isapprox(σ(rp_res[NeXLMatrixCorrection.RLabel(matname,inner(cxr))]), 0.00146, atol=0.001)
@test isapprox(σ(rp_res[NeXLMatrixCorrection.ϕ0Label(matname, inner(cxr))]), 0.00643, atol=0.0001)

maintain = MaintainLabels([NeXLMatrixCorrection.ZbarbLabel, NeXLMatrixCorrection.ϕ0Label, NeXLMatrixCorrection.E0keVLabel], rp_res)

frbar = NeXLMatrixCorrection.StepFRBar(matname, inner(cxr)) | maintain
frbar_res = frbar(rp_res)
frbar_model = frbar ∘ rp_model
frbar_mcres = mcpropagate(frbar_model, k240, 1000, parallel=false, rng=rgen)

#println("Analytical Result")
#print(frbar_res)
#println("MC Result")
#print(frbar_mcres)

@test isapprox(value(frbar_res[NeXLMatrixCorrection.FLabel(matname,inner(cxr))]), xpp.F, atol=0.01*xpp.F)
@test isapprox(value(frbar_res[NeXLMatrixCorrection.RbarLabel(matname,inner(cxr))]), 3.31e-4, atol=0.01e-4)
@test isapprox(σ(frbar_res[NeXLMatrixCorrection.FLabel(matname,inner(cxr))]), 2.56e-5, atol=0.1e-5)
@test isapprox(σ(frbar_res[NeXLMatrixCorrection.RbarLabel(matname, inner(cxr))]), 5.36e-6, atol=0.1e-6)

maintain = MaintainLabels([NeXLMatrixCorrection.ϕ0Label, NeXLMatrixCorrection.RbarLabel, NeXLMatrixCorrection.FLabel], frbar_res)

pb = NeXLMatrixCorrection.StepPb(matname, inner(cxr)) | maintain
pb_res = pb(frbar_res)
pb_model = pb ∘ frbar_model
pb_mcres = mcpropagate(pb_model, k240, 1000, parallel=false, rng=rgen)

@test isapprox(value(pb_res[NeXLMatrixCorrection.PLabel(matname,inner(cxr))]), 1.42e4, atol=0.01e4)
@test isapprox(value(pb_res[NeXLMatrixCorrection.bLabel(matname,inner(cxr))]), 7.78e3, atol=0.01e3)
@test isapprox(σ(pb_res[NeXLMatrixCorrection.PLabel(matname,inner(cxr))]), 2.85e2, atol=0.4e2)
@test isapprox(σ(pb_res[NeXLMatrixCorrection.bLabel(matname,inner(cxr))]), 1.25e2, atol=0.2e2)

maintain = MaintainLabels([NeXLMatrixCorrection.ϕ0Label, NeXLMatrixCorrection.FLabel, NeXLMatrixCorrection.PLabel, NeXLMatrixCorrection.bLabel], pb_res)

aϵ = NeXLMatrixCorrection.Stepaϵ(matname, inner(cxr)) | maintain
aϵ_res = aϵ(pb_res)
aϵ_model = aϵ ∘ pb_model
aϵ_mcres = mcpropagate(aϵ_model, k240, 1000, parallel=false, rng=rgen)

@test isapprox(value(aϵ_res[NeXLMatrixCorrection.aLabel(matname,inner(cxr))]), 6.75e3, atol=0.01e3)
@test isapprox(value(aϵ_res[NeXLMatrixCorrection.ϵLabel(matname,inner(cxr))]), -0.133, atol=0.001)
@test isapprox(σ(aϵ_res[NeXLMatrixCorrection.aLabel(matname,inner(cxr))]), 100, atol=10)
@test isapprox(σ(aϵ_res[NeXLMatrixCorrection.ϵLabel(matname,inner(cxr))]), 0.0013, atol=0.0001)

maintain = MaintainLabels( [NeXLMatrixCorrection.bLabel, NeXLMatrixCorrection.ϕ0Label, NeXLMatrixCorrection.ϵLabel ], aϵ_res)

AB = NeXLMatrixCorrection.StepAB(matname, inner(cxr)) | maintain
AB_res = AB(aϵ_res)
AB_model = AB ∘ aϵ_model
AB_mcres = mcpropagate(AB_model, k240, 1000, parallel=false, rng=rgen)

@test isapprox(value(AB_res[NeXLMatrixCorrection.ALabel(matname,inner(cxr))]), 393, atol=1)
@test isapprox(value(AB_res[NeXLMatrixCorrection.BLabel(matname,inner(cxr))]), -3.79e5, atol=0.01e5)
@test isapprox(σ(AB_res[NeXLMatrixCorrection.ALabel(matname,inner(cxr))]), 8.46, atol=0.01)
@test isapprox(σ(AB_res[NeXLMatrixCorrection.BLabel(matname,inner(cxr))]), 4.0e3, atol=4.0e3)

χdata =  uvs(
    NeXLMatrixCorrection.μoρLabel(matname,inner(cxr))=>uv(5813.,577.), #
    NeXLMatrixCorrection.θLabel(matname)=>uv(deg2rad(40.0), deg2rad(0.1)), #
    NeXLMatrixCorrection.dzLabel(matname)=>uv(0.0, 1.0e-6) # 10 nm
)

ABχ_res = cat(AB_res, χdata)

χFr = NeXLMatrixCorrection.StepχFr(matname, inner(cxr))
χFr_res = χFr(ABχ_res)
χFr_model = χFr ∘ AB_model
χFr_mcres = mcpropagate(χFr, cat(AB_res,χdata), 1000, parallel=false, rng=rgen)


#end
