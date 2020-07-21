using NeXLCore
using NeXLMatrixCorrection
using Random
using Test

@testset "XPPU" begin
    rgen = MersenneTwister(0xBADF00D)

    unknown, standard, cxr, e0, θtoa = "K240", "SiO2", n"O K-L3", 20.0e3, deg2rad(40.0)
    unk_mat = material(
        unknown, #
        n"O" => 0.340023,
        n"Mg" => 0.030154,
        n"Si" => 0.186986,
        n"Ti" => 0.059950, #
        n"Zn" => 0.040168,
        n"Zr" => 0.074030,
        n"Ba" => 0.268689,
    )
    std_mat = parse(Material, "SiO2")

    # For comparison with direct calculation...
    zaf = zafcorrection(
        XPP,
        NullFluorescence,
        Coating,
        unk_mat,
        std_mat,
        inner(cxr),
        e0,
        unkCoating = carboncoating(10.0),
        stdCoating = carboncoating(11.0),
    )

    m = NeXLMatrixCorrection.m(XPP, inner(cxr))

    input_uvs = uvs(
        (MassFractionLabel(unknown, elm) => uv(unk_mat[elm], 0.01 * unk_mat[elm]) for elm in keys(unk_mat))...,
        (AtomicWeightLabel(unknown, elm) => uv(a(elm, unk_mat), 0.001 * a(elm, unk_mat)) for elm in keys(unk_mat))...,
        (JzLabel(elm) => Ju(elm) for elm in keys(unk_mat))...,
        NeXLMatrixCorrection.E0Label(unknown) => uv(e0, 0.1e3),
        NeXLMatrixCorrection.mLabel(inner(cxr)) => uv(m, 0.01 * m),
    )

    maintain = MaintainInputs([NeXLMatrixCorrection.mLabel], input_uvs)
    mjz = StepMJZbarb(unknown, [n"O", n"Mg", n"Si", n"Ti", n"Zn", n"Zr", n"Ba"]) | maintain

    mjz_res = mjz(input_uvs)
    mjz_mcres = mcpropagate(mjz, input_uvs, 10000, parallel = false, rng = rgen)

    @test isapprox(value(mjz_res[NeXLMatrixCorrection.BigMLabel(unknown)]), 0.466, atol = 0.001)
    @test isapprox(value(mjz_res[NeXLMatrixCorrection.ZbarbLabel(unknown)]), 22.5, atol = 0.1)
    @test isapprox(value(mjz_res[NeXLMatrixCorrection.E0keVLabel(unknown)]), 20.0, atol = 0.1)
    @test isapprox(value(mjz_res[NeXLMatrixCorrection.JLabel(unknown)]), 0.217, atol = 0.001)
    @test isapprox(σ(mjz_res[NeXLMatrixCorrection.BigMLabel(unknown)]), 0.0023, atol = 0.0001)
    @test isapprox(σ(mjz_res[NeXLMatrixCorrection.ZbarbLabel(unknown)]), 0.25, atol = 0.03)
    @test isapprox(σ(mjz_res[NeXLMatrixCorrection.E0keVLabel(unknown)]), 0.1, atol = 0.01)
    @test isapprox(σ(mjz_res[NeXLMatrixCorrection.JLabel(unknown)]), 0.00198, atol = 0.0001)

    # println("Analytical Result")
    # print(mjz_res)
    # println("MC Result")
    # print(mjz_mcres)

    dpt = StepDPT(unknown, inner(cxr)) | AllInputs()

    dpt_res = dpt(mjz_res)
    dpt_mcres = mcpropagate(dpt ∘ mjz, input_uvs, 1000, parallel = false, rng = rgen)

    @test isapprox(value(dpt_res[NeXLMatrixCorrection.DkLabel(unknown, inner(cxr), 1)]), 6.6e-6, atol = 0.01e-6)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.DkLabel(unknown, inner(cxr), 2)]), 1.49e-05, atol = 0.01e-5)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.DkLabel(unknown, inner(cxr), 3)]), 1.02e-5, atol = 0.01e-5)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.PkLabel(unknown, inner(cxr), 1)]), 0.780, atol = 0.001)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.PkLabel(unknown, inner(cxr), 2)]), 0.100, atol = 0.001)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.PkLabel(unknown, inner(cxr), 3)]), -0.446, atol = 0.001)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.TkLabel(unknown, inner(cxr), 1)]), 0.911, atol = 0.001)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.TkLabel(unknown, inner(cxr), 2)]), 0.231, atol = 0.001)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.TkLabel(unknown, inner(cxr), 3)]), -0.315, atol = 0.001)

    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.DkLabel(unknown, inner(cxr), 1)]), 0.0, atol = 0.001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.DkLabel(unknown, inner(cxr), 2)]), 4.3e-09, atol = 0.1e-9)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.DkLabel(unknown, inner(cxr), 3)]), 9.3e-08, atol = 0.1e-8)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.PkLabel(unknown, inner(cxr), 1)]), 0.0, atol = 0.001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.PkLabel(unknown, inner(cxr), 2)]), 0.0, atol = 0.001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.PkLabel(unknown, inner(cxr), 3)]), 0.00049, atol = 0.00001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.TkLabel(unknown, inner(cxr), 1)]), 0.009, atol = 0.001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.TkLabel(unknown, inner(cxr), 2)]), 0.009, atol = 0.001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.TkLabel(unknown, inner(cxr), 3)]), 0.0088, atol = 0.0001)


    @test isapprox(value(dpt_res[NeXLMatrixCorrection.BigMLabel(unknown)]), 0.466, atol = 0.001)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.ZbarbLabel(unknown)]), 22.5, atol = 0.1)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.E0keVLabel(unknown)]), 20.0, atol = 0.1)
    @test isapprox(value(dpt_res[NeXLMatrixCorrection.JLabel(unknown)]), 0.217, atol = 0.001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.BigMLabel(unknown)]), 0.0023, atol = 0.0001)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.ZbarbLabel(unknown)]), 0.25, atol = 0.03)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.E0keVLabel(unknown)]), 0.1, atol = 0.01)
    @test isapprox(σ(dpt_res[NeXLMatrixCorrection.JLabel(unknown)]), 0.00198, atol = 0.0001)

    # println("Analytical Result")
    # print(dpt_res)
    # println("MC Result")
    # print(dpt_mcres)

    retain = MaintainInputs([NeXLMatrixCorrection.E0keVLabel, NeXLMatrixCorrection.ZbarbLabel], dpt_res)

    qla = StepQlaOoS(unknown, inner(cxr)) | retain
    qla_res = qla(dpt_res)
    qla_model = qla ∘ dpt ∘ mjz
    qla_mcres = mcpropagate(qla_model, input_uvs, 10000, parallel = false, rng = rgen)

    #println("Analytical Result")
    #print(qla_res)
    #println("MC Result")
    #print(qla_mcres)

    @test isapprox(value(qla_res[NeXLMatrixCorrection.QlaLabel(unknown, inner(cxr))]), 0.548, atol = 0.001)
    @test isapprox(value(qla_res[NeXLMatrixCorrection.OoSLabel(unknown, inner(cxr))]), 0.00108, atol = 0.00001)
    @test isapprox(value(qla_res[NeXLMatrixCorrection.ηLabel(unknown)]), 0.253, atol = 0.001)
    @test isapprox(value(qla_res[NeXLMatrixCorrection.JU0Label(unknown)]), 99.8, atol = 0.1)
    @test isapprox(value(qla_res[NeXLMatrixCorrection.WbarLabel(unknown)]), 0.665, atol = 0.001)
    @test isapprox(value(qla_res[NeXLMatrixCorrection.qLabel(unknown)]), 0.987, atol = 0.001)

    @test isapprox(σ(qla_res[NeXLMatrixCorrection.QlaLabel(unknown, inner(cxr))]), 0.0173, atol = 0.001)
    @test isapprox(σ(qla_res[NeXLMatrixCorrection.OoSLabel(unknown, inner(cxr))]), 2.85e-5, atol = 0.1e-5)
    @test isapprox(σ(qla_res[NeXLMatrixCorrection.ηLabel(unknown)]), 0.00218, atol = 0.0001)
    @test isapprox(σ(qla_res[NeXLMatrixCorrection.JU0Label(unknown)]), 0.682, atol = 0.001)
    @test isapprox(σ(qla_res[NeXLMatrixCorrection.WbarLabel(unknown)]), 0.000667, atol = 0.00001)
    @test isapprox(σ(qla_res[NeXLMatrixCorrection.qLabel(unknown)]), 0.00595, atol = 0.0001)

    maintain = MaintainInputs(
        [
            NeXLMatrixCorrection.E0keVLabel,
            NeXLMatrixCorrection.ZbarbLabel,
            NeXLMatrixCorrection.OoSLabel,
            NeXLMatrixCorrection.QlaLabel,
        ],
        qla_res,
    )

    rp = NeXLMatrixCorrection.StepRPhi0(unknown, inner(cxr)) | maintain
    rp_res = rp(qla_res)
    rp_model = rp ∘ qla_model
    rp_mcres = mcpropagate(rp_model, input_uvs, 1000, parallel = false, rng = rgen)

    #println("Analytical Result")
    #print(rp_res)
    #println("MC Result")
    #print(rp_mcres)

    @test isapprox(value(rp_res[NeXLMatrixCorrection.RLabel(unknown, inner(cxr))]), 0.852, atol = 0.001)
    @test isapprox(value(rp_res[NeXLMatrixCorrection.ϕ0Label(unknown, inner(cxr))]), zaf[1].za.ϕ0, atol = 1.0e-7)

    @test isapprox(σ(rp_res[NeXLMatrixCorrection.RLabel(unknown, inner(cxr))]), 0.00146, atol = 0.001)
    @test isapprox(σ(rp_res[NeXLMatrixCorrection.ϕ0Label(unknown, inner(cxr))]), 0.00643, atol = 0.0001)

    maintain = MaintainInputs(
        [NeXLMatrixCorrection.ZbarbLabel, NeXLMatrixCorrection.ϕ0Label, NeXLMatrixCorrection.E0keVLabel],
        rp_res,
    )

    frbar = NeXLMatrixCorrection.StepFRBar(unknown, inner(cxr)) | maintain
    frbar_res = frbar(rp_res)
    frbar_model = frbar ∘ rp_model
    frbar_mcres = mcpropagate(frbar_model, input_uvs, 1000, parallel = false, rng = rgen)

    #println("Analytical Result")
    #print(frbar_res)
    #println("MC Result")
    #print(frbar_mcres)

    @test isapprox(value(frbar_res[NeXLMatrixCorrection.FLabel(unknown, inner(cxr))]), 0.00168, atol = 0.00001)
    @test isapprox(value(frbar_res[NeXLMatrixCorrection.RbarLabel(unknown, inner(cxr))]), 3.31e-4, atol = 0.01e-4)
    @test isapprox(σ(frbar_res[NeXLMatrixCorrection.FLabel(unknown, inner(cxr))]), 2.56e-5, atol = 0.1e-5)
    @test isapprox(σ(frbar_res[NeXLMatrixCorrection.RbarLabel(unknown, inner(cxr))]), 5.36e-6, atol = 0.1e-6)

    @test isapprox(value(frbar_res[NeXLMatrixCorrection.FLabel(unknown, inner(cxr))]), zaf[1].za.F, rtol = 1.0e-7)


    maintain = MaintainInputs(
        [NeXLMatrixCorrection.ϕ0Label, NeXLMatrixCorrection.RbarLabel, NeXLMatrixCorrection.FLabel],
        frbar_res,
    )
    pb = NeXLMatrixCorrection.StepPb(unknown, inner(cxr)) | maintain
    pb_res = pb(frbar_res)
    pb_model = pb ∘ frbar_model
    pb_mcres = mcpropagate(pb_model, input_uvs, 1000, parallel = false, rng = rgen)

    @test isapprox(value(pb_res[NeXLMatrixCorrection.PLabel(unknown, inner(cxr))]), 1.42e4, atol = 0.01e4)
    @test isapprox(value(pb_res[NeXLMatrixCorrection.bLabel(unknown, inner(cxr))]), 7.78e3, atol = 0.01e3)
    @test isapprox(σ(pb_res[NeXLMatrixCorrection.PLabel(unknown, inner(cxr))]), 2.85e2, atol = 0.4e2)
    @test isapprox(σ(pb_res[NeXLMatrixCorrection.bLabel(unknown, inner(cxr))]), 1.25e2, atol = 0.2e2)

    maintain = MaintainInputs(
        [
            NeXLMatrixCorrection.ϕ0Label,
            NeXLMatrixCorrection.FLabel,
            NeXLMatrixCorrection.PLabel,
            NeXLMatrixCorrection.bLabel,
        ],
        pb_res,
    )
    aϵ = NeXLMatrixCorrection.Stepaϵ(unknown, inner(cxr)) | maintain
    aϵ_res = aϵ(pb_res)
    aϵ_model = aϵ ∘ pb_model
    aϵ_mcres = mcpropagate(aϵ_model, input_uvs, 1000, parallel = false, rng = rgen)

    @test isapprox(value(aϵ_res[NeXLMatrixCorrection.aLabel(unknown, inner(cxr))]), 6.75e3, atol = 0.01e3)
    @test isapprox(value(aϵ_res[NeXLMatrixCorrection.aLabel(unknown, inner(cxr))]), zaf[1].za.a, rtol = 1.0e-7)
    @test isapprox(value(aϵ_res[NeXLMatrixCorrection.ϵLabel(unknown, inner(cxr))]), -0.133, atol = 0.001)
    @test isapprox(σ(aϵ_res[NeXLMatrixCorrection.aLabel(unknown, inner(cxr))]), 100, atol = 10)
    @test isapprox(σ(aϵ_res[NeXLMatrixCorrection.ϵLabel(unknown, inner(cxr))]), 0.0013, atol = 0.0001)

    maintain = MaintainInputs(
        [
            NeXLMatrixCorrection.bLabel,
            NeXLMatrixCorrection.ϕ0Label,
            NeXLMatrixCorrection.ϵLabel,
            NeXLMatrixCorrection.FLabel,
        ],
        aϵ_res,
    )
    AB = NeXLMatrixCorrection.StepAB(unknown, inner(cxr)) | maintain
    AB_res = AB(aϵ_res)
    AB_model = AB ∘ aϵ_model
    AB_mcres = mcpropagate(AB_model, input_uvs, 1000, parallel = false, rng = rgen)

    @test isapprox(value(AB_res[NeXLMatrixCorrection.ALabel(unknown, inner(cxr))]), 393, atol = 1)
    @test isapprox(value(AB_res[NeXLMatrixCorrection.BLabel(unknown, inner(cxr))]), -3.79e5, atol = 0.01e5)
    @test isapprox(σ(AB_res[NeXLMatrixCorrection.ALabel(unknown, inner(cxr))]), 8.46, atol = 0.01)
    @test isapprox(σ(AB_res[NeXLMatrixCorrection.BLabel(unknown, inner(cxr))]), 4.0e3, atol = 4.0e3)

    @test isapprox(value(AB_res[NeXLMatrixCorrection.ALabel(unknown, inner(cxr))]), zaf[1].za.A, rtol = 1.0e-7)
    @test isapprox(value(AB_res[NeXLMatrixCorrection.BLabel(unknown, inner(cxr))]), zaf[1].za.B, rtol = 1.0e-7)


    χdata = uvs(
        NeXLMatrixCorrection.μoρLabel(unknown, cxr) => uv(5813.0, 577.0), #
        NeXLMatrixCorrection.θLabel(unknown) => uv(θtoa, deg2rad(0.1)), #
        NeXLMatrixCorrection.dzLabel(unknown) => uv(0.0, 1.0e-6), # 10 nm
    )

    ABχ_res = cat(AB_res, χdata)

    maintain = MaintainInputs([NeXLMatrixCorrection.θLabel, NeXLMatrixCorrection.FLabel], ABχ_res)
    χFr = NeXLMatrixCorrection.StepχFr(unknown, cxr) | maintain
    χFr_res = χFr(ABχ_res)
    χFr_model = χFr ∘ AB_model
    # χFr_mcres = mcpropagate(χFr, ABχ_res, 1000, parallel=false, rng=rgen)

    @test isapprox(value(χFr_res[NeXLMatrixCorrection.χLabel(unknown, cxr)]), 9.04e3, atol = 0.01e3)
    @test isapprox(value(χFr_res[NeXLMatrixCorrection.FrLabel(unknown, cxr)]), 2.85e-4, atol = 0.1e-4)
    @test isapprox(σ(χFr_res[NeXLMatrixCorrection.χLabel(unknown, cxr)]), 898, atol = 10)
    @test isapprox(σ(χFr_res[NeXLMatrixCorrection.FrLabel(unknown, cxr)]), 2.6e-5, atol = 0.1e-5)

    coatU = "10 nm C"
    coatingData = uvs( #
        NeXLMatrixCorrection.tcLabel(coatU) => uv(1.0e-6, 0.5e-6),
        NeXLMatrixCorrection.μoρLabel(coatU, cxr) => uv(11705.0, 2185.0),  # n"O K-L3" in n"C"
    )

    Frc_input = cat(χFr_res, coatingData)

    maintain = MaintainInputs([NeXLMatrixCorrection.θLabel, NeXLMatrixCorrection.FLabel], Frc_input)
    Frc = NeXLMatrixCorrection.StepFrc(unknown, coatU, cxr) | maintain
    Frcu_res = Frc(Frc_input)
    Frc_model = Frc
    Frc_mcres = mcpropagate(Frc, Frc_input, 1000, parallel = false, rng = rgen)

    @test isapprox(value(Frcu_res[NeXLMatrixCorrection.FrcLabel(unknown, coatU, cxr)]), 2.8e-4, atol = 0.1e-4)
    @test isapprox(σ(Frcu_res[NeXLMatrixCorrection.FrcLabel(unknown, coatU, cxr)]), 2.5e-5, atol = 0.1e-5)

    # Now calculate the standard

    input_uvs = uvs(
        (MassFractionLabel(standard, elm) => uv(std_mat[elm], 0.01 * std_mat[elm]) for elm in keys(std_mat))...,
        (AtomicWeightLabel(standard, elm) => uv(a(elm, std_mat), 0.001 * a(elm, std_mat)) for elm in keys(std_mat))...,
        (JzLabel(elm) => Ju(elm) for elm in keys(std_mat))...,
        NeXLMatrixCorrection.E0Label(standard) => uv(e0, 0.1e3),
        NeXLMatrixCorrection.mLabel(inner(cxr)) => uv(m, 0.01 * m),
    )

    maintain = MaintainInputs([NeXLMatrixCorrection.mLabel], input_uvs)
    mjz = StepMJZbarb(standard, [n"O", n"Si"]) | maintain

    mjz_res = mjz(input_uvs)
    mjz_mcres = mcpropagate(mjz, input_uvs, 10000, parallel = false, rng = rgen)

    # println("Analytical Result")
    # print(mjz_res)
    # println("MC Result")
    # print(mjz_mcres)

    dpt = StepDPT(standard, inner(cxr)) | AllInputs()

    dpt_res = dpt(mjz_res)
    dpt_mcres = mcpropagate(dpt ∘ mjz, input_uvs, 1000, parallel = false, rng = rgen)

    # println("Analytical Result")
    # print(dpt_res)
    # println("MC Result")
    # print(dpt_mcres)

    retain = MaintainInputs([NeXLMatrixCorrection.E0keVLabel, NeXLMatrixCorrection.ZbarbLabel], dpt_res)
    qla = StepQlaOoS(standard, inner(cxr)) | retain
    qla_res = qla(dpt_res)
    qla_model = qla ∘ dpt ∘ mjz
    qla_mcres = mcpropagate(qla_model, input_uvs, 10000, parallel = false, rng = rgen)

    #println("Analytical Result")
    #print(qla_res)
    #println("MC Result")
    #print(qla_mcres)

    maintain = MaintainInputs(
        [
            NeXLMatrixCorrection.E0keVLabel,
            NeXLMatrixCorrection.ZbarbLabel,
            NeXLMatrixCorrection.OoSLabel,
            NeXLMatrixCorrection.QlaLabel,
        ],
        qla_res,
    )
    rp = NeXLMatrixCorrection.StepRPhi0(standard, inner(cxr)) | maintain
    rp_res = rp(qla_res)
    rp_model = rp ∘ qla_model
    rp_mcres = mcpropagate(rp_model, input_uvs, 1000, parallel = false, rng = rgen)

    #println("Analytical Result")
    #print(rp_res)
    #println("MC Result")
    #print(rp_mcres)

    maintain = MaintainInputs(
        [NeXLMatrixCorrection.ZbarbLabel, NeXLMatrixCorrection.ϕ0Label, NeXLMatrixCorrection.E0keVLabel],
        rp_res,
    )
    frbar = NeXLMatrixCorrection.StepFRBar(standard, inner(cxr)) | maintain
    frbar_res = frbar(rp_res)
    frbar_model = frbar ∘ rp_model
    frbar_mcres = mcpropagate(frbar_model, input_uvs, 1000, parallel = false, rng = rgen)

    #println("Analytical Result")
    #print(frbar_res)
    #println("MC Result")
    #print(frbar_mcres)

    maintain = MaintainInputs(
        [NeXLMatrixCorrection.ϕ0Label, NeXLMatrixCorrection.RbarLabel, NeXLMatrixCorrection.FLabel],
        frbar_res,
    )
    pb = NeXLMatrixCorrection.StepPb(standard, inner(cxr)) | maintain
    pb_res = pb(frbar_res)
    pb_model = pb ∘ frbar_model
    pb_mcres = mcpropagate(pb_model, input_uvs, 1000, parallel = false, rng = rgen)

    maintain = MaintainInputs(
        [
            NeXLMatrixCorrection.ϕ0Label,
            NeXLMatrixCorrection.FLabel,
            NeXLMatrixCorrection.PLabel,
            NeXLMatrixCorrection.bLabel,
        ],
        pb_res,
    )
    aϵ = NeXLMatrixCorrection.Stepaϵ(standard, inner(cxr)) | maintain
    aϵ_res = aϵ(pb_res)
    aϵ_model = aϵ ∘ pb_model
    aϵ_mcres = mcpropagate(aϵ_model, input_uvs, 1000, parallel = false, rng = rgen)

    maintain = MaintainInputs(
        [
            NeXLMatrixCorrection.bLabel,
            NeXLMatrixCorrection.ϕ0Label,
            NeXLMatrixCorrection.ϵLabel,
            NeXLMatrixCorrection.FLabel,
        ],
        aϵ_res,
    )
    AB = NeXLMatrixCorrection.StepAB(standard, inner(cxr)) | maintain
    AB_res = AB(aϵ_res)
    AB_model = AB ∘ aϵ_model
    AB_mcres = mcpropagate(AB_model, input_uvs, 1000, parallel = false, rng = rgen)

    χdata = uvs(
        NeXLMatrixCorrection.μoρLabel(standard, cxr) => uv(4123.0, 668.0), #
        NeXLMatrixCorrection.θLabel(standard) => uv(θtoa, deg2rad(0.1)), #
        NeXLMatrixCorrection.dzLabel(standard) => uv(0.0, 1.0e-6), # 10 nm
    )

    ABχ_res = cat(AB_res, χdata)

    maintain = MaintainInputs([NeXLMatrixCorrection.θLabel, NeXLMatrixCorrection.FLabel], ABχ_res)
    χFr = NeXLMatrixCorrection.StepχFr(standard, cxr) | maintain
    χFr_res = χFr(ABχ_res)
    χFr_model = χFr ∘ AB_model
    input = cat(AB_res, χdata)
    display(input)
    println()
    χFr_mcres = mcpropagate(χFr, input, 1000, parallel = false, rng = rgen)

    coatS = "11 nm C"
    coatingData = uvs( #
        NeXLMatrixCorrection.tcLabel(coatS) => uv(1.1e-6, 0.5e-6),
        NeXLMatrixCorrection.μoρLabel(coatS, cxr) => uv(mac(pure(n"C"), cxr), 2185.0),  # n"O K-L3" in n"C"
    )
    Frc_input = cat(χFr_res, coatingData)

    maintain = MaintainInputs([NeXLMatrixCorrection.θLabel, NeXLMatrixCorrection.FLabel], Frc_input)
    Frc = NeXLMatrixCorrection.StepFrc(standard, coatS, cxr) | maintain
    Frcs_res = Frc(Frc_input)
    Frc_model = Frc
    Frc_mcres = mcpropagate(Frc, Frc_input, 1000, parallel = false, rng = rgen)

    za_inputs = cat(Frcu_res, Frcs_res)

    za = NeXLMatrixCorrection.StepZA(unknown, standard, cxr, coatU, coatS)
    za_res = za(za_inputs)
    za_mcres = mcpropagate(za, za_inputs, 1000, parallel = false, rng = rgen)

    # Compare to DTSA-II and MC
    @test isapprox(value(za_res[NeXLMatrixCorrection.ZLabel(unknown, standard, inner(cxr))]), 1.1267, atol = 0.001)
    @test isapprox(
        value(za_res[NeXLMatrixCorrection.AbsLabel(unknown, standard, cxr, coatU, coatS)]),
        0.762,
        atol = 0.01,
    )
    @test isapprox(σ(za_res[NeXLMatrixCorrection.ZLabel(unknown, standard, inner(cxr))]), 0.024, atol = 0.001)
    @test isapprox(σ(za_res[NeXLMatrixCorrection.AbsLabel(unknown, standard, cxr, coatU, coatS)]), 0.14, atol = 0.01)

    zaf_za = ZAFc(zaf..., cxr, θtoa, θtoa)
    zaf_a = ZAFc(zaf..., cxr, θtoa, θtoa) / Z(zaf...)

    @test isapprox(value(za_res[NeXLMatrixCorrection.ZLabel(unknown, standard, inner(cxr))]), Z(zaf...), atol = 1e-5)
    @test_broken isapprox(
        value(za_res[NeXLMatrixCorrection.AbsLabel(unknown, standard, cxr, coatU, coatS)]),
        zaf_a,
        atol = 1e-5,
    )
    @test_broken isapprox(
        value(za_res[NeXLMatrixCorrection.ZALabel(unknown, standard, cxr, coatU, coatS)]),
        zaf_za,
        atol = 1e-5,
    )
end
