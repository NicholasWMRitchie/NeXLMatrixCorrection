using NeXLUncertainties

asserting() = false

# Removeable asserts
macro Assert(test)
    esc(:(
        if $(@__MODULE__).asserting()
            @assert($test)
        end
    ))
end

# Note:: XPP is detailed in the "Green Book" Electron Probe Quantiication in terms using keV for energy.  Since
#   NeXL uses eV by default, I convert all inputs to keV in the first step and use keV thoughout this calculation.


struct JzLabel <: Label  # in eV
    element::Element
end

Base.show(io::IO, jz::JzLabel) = print(io, "Jz[", jz.element.symbol, "]")

struct StepMJZbarb <: MeasurementModel
    material::String
    elements::Vector{Element}
end

Base.show(io::IO, mjz::StepMJZbarb) = print(io, "MJZbarb[", mjz.material, "]")

struct BigMLabel <: Label
    material::String
end
Base.show(io::IO, m::BigMLabel) = print(io, "M[", m.material, "]")

struct JLabel <: Label # in keV
    material::String
end
Base.show(io::IO, j::JLabel) = print(io, "J[", j.material, "]")

struct E0Label <: Label
    material::String
end
Base.show(io::IO, e0::E0Label) = print(io, "E0[", e0.material, "]")

struct E0keVLabel <: Label
    material::String
end
Base.show(io::IO, e0::E0keVLabel) = print(io, "E0(keV)[", e0.material, "]")

struct ZbarbLabel <: Label
    material::String
end
Base.show(io::IO, zb::ZbarbLabel) = print(io, "Zbarb[", zb.material, "]")

function Ju(elm::Element, f = 0.01) #C1
    j = NeXLMatrixCorrection.J(Zeller1973, elm) # from xpp.jl (in eV)
    return uv(j, f * j)
end

function NeXLUncertainties.compute(mjz::StepMJZbarb, inputs::LabeledValues, withJac::Bool)::MMResult
    # Build the labels once...
    mfls = map(elm -> MassFractionLabel(mjz.material, elm), mjz.elements)
    awls = map(elm -> AtomicWeightLabel(mjz.material, elm), mjz.elements)
    jzs = map(elm -> JzLabel(elm), mjz.elements)
    e0l = E0Label(mjz.material)
    # Unpack all the variables
    c, a, = map(mfl -> inputs[mfl], mfls), map(awl -> inputs[awl], awls)
    z = map(elm -> convert(Float64, elm.number), mjz.elements)
    keV = 0.001
    j = map(jz -> keV * inputs[jz], jzs) # in keV
    e0 = inputs[e0l]
    # Outputs...
    outputs = [BigMLabel(mjz.material), JLabel(mjz.material), ZbarbLabel(mjz.material), E0keVLabel(mjz.material)]
    M = mapreduce(i -> c[i] * z[i] / a[i], +, eachindex(z))
    J = exp(sum((c[i] * z[i] / a[i]) * log(j[i]) for i in eachindex(z)) / M) # in keV
    Zb = sum(c[i] * sqrt(z[i]) for i in eachindex(z))^2
    values = [M, J, Zb, keV * e0]
    jacob = withJac ? zeros(Float64, length(values), length(inputs)) : missing
    if withJac
        for i in eachindex(z)
            mfl, awl, jz = mfls[i], awls[i], jzs[i]
            # dM/d? (index = 1 for M)
            jacob[1, indexin(mfl, inputs)] = z[i] / a[i]
            jacob[1, indexin(awl, inputs)] = -c[i] * z[i] / (a[i]^2)
            # dJ/d? (index = 2 for J)
            kk = (j[i] * z[i]) / (M * a[i])
            jacob[2, indexin(jz, inputs)] = keV * kk * (c[i] / j[i])
            jacob[2, indexin(mfl, inputs)] = kk * (log(j[i]) - log(J))
            jacob[2, indexin(awl, inputs)] = kk * (c[i] / a[i]) * (log(J) - log(j[i]))
            # dZ/d? (index = 3 for Zbarb)
            jacob[3, indexin(mfl, inputs)] = 2.0 * sqrt(z[i] * Zb)
            # dE0kevdE0 (index = 4 for E0keV)
            jacob[4, indexin(e0l, inputs)] = keV
        end
    end
    return (LabeledValues(outputs, values), jacob)
end

struct StepDPT <: MeasurementModel
    material::String
    shell::AtomicSubShell
end

struct mLabel <: Label
    shell::AtomicSubShell
end
Base.show(io::IO, m::mLabel) = print(io, "m[", m.shell, "]")

struct DkLabel <: Label
    material::String
    shell::AtomicSubShell
    k::Int
end
Base.show(io::IO, d::DkLabel) = print(io, "D[", d.material, ",", d.shell, ",", d.k, "]")

struct PkLabel <: Label
    material::String
    shell::AtomicSubShell
    k::Int
end
Base.show(io::IO, p::PkLabel) = print(io, "P[", p.material, ",", p.shell, ",", p.k, "]")

struct TkLabel <: Label
    material::String
    shell::AtomicSubShell
    k::Int
end
Base.show(io::IO, t::TkLabel) = print(io, "T[", t.material, ",", t.shell, ",", t.k, "]")

function NeXLUncertainties.compute(dpt::StepDPT, inputs::LabeledValues, withJac::Bool)::MMResult
    # inputs
    Jl, ml = JLabel(dpt.material), mLabel(dpt.shell)
    J, m = inputs[Jl], inputs[ml]
    # outputs
    Dls = [DkLabel(dpt.material, dpt.shell, k) for k = 1:3]
    Pls = [PkLabel(dpt.material, dpt.shell, k) for k = 1:3]
    Tls = [TkLabel(dpt.material, dpt.shell, k) for k = 1:3]

    D = (6.6e-6, 1.12e-5 * (1.35 - 0.45 * J^2), 2.2e-6 / J)
    P = (0.78, 0.1, 0.25J - 0.5)
    T = (1.0 - m) .+ P

    vals = LabeledValues([Dls..., Pls..., Tls...], [D..., P..., T...])
    jacob = withJac ? zeros(Float64, length(vals), length(inputs)) : missing
    if withJac
        Jli, mli = indexin(Jl, inputs), indexin(ml, inputs)
        # dD/dJ = ( 0.0, 2.0*1.12e-5*-0.45*J), -2.2e-6/(J^2) ))
        @Assert indexin(Dls[2], vals) == 2
        jacob[2, Jli] = -1.008e-5 * J   # indexin(Dls[2],vals)==2
        @Assert indexin(Dls[3], vals) == 3
        jacob[3, Jli] = (-2.2e-6 / (J^2))# indexin(Dls[3],vals)==3
        # dP/dJ = dP/dJ = ( 0.0, 0.0, 0.25 )
        @Assert indexin(Pls[3], vals) == 3 + 3
        jacob[3+3, Jli] = 0.25           # indexin(Pls[3],vals) == 3+3
        @Assert indexin(Tls[3], vals) == 6 + 3
        jacob[6+3, Jli] = 0.25           # indexin(Tls[3],vals) == 6+3
        # dT/dm = ( -1.0, -1.0, -1.0 )
        for k = 1:3
            @Assert indexin(Tls[k], vals) == 6 + k
            jacob[6+k, mli] = -1.0    # indexin(Tls[k],vals)==6+k
        end
    end
    return (vals, jacob)
end

struct StepQlaOoS <: MeasurementModel
    material::String
    shell::AtomicSubShell
end

struct QlaLabel <: Label
    material::String
    shell::AtomicSubShell
end
Base.show(io::IO, qla::QlaLabel) = print(io, "Qla[", qla.material, ",", qla.shell, "]")

struct ηLabel <: Label
    material::String
end
Base.show(io::IO, m::ηLabel) = print(io, "η[", m.material, "]")

struct JU0Label <: Label
    material::String
end
Base.show(io::IO, m::JU0Label) = print(io, "J(U0)[", m.material, "]")

struct WbarLabel <: Label
    material::String
end
Base.show(io::IO, m::WbarLabel) = print(io, "Wbar[", m.material, "]")

struct qLabel <: Label
    material::String
end
Base.show(io::IO, m::qLabel) = print(io, "q[", m.material, "]")

struct OoSLabel <: Label
    material::String
    shell::AtomicSubShell
end
Base.show(io::IO, d::OoSLabel) = print(io, "1/S[", d.material, ",", d.shell, "]")

function NeXLUncertainties.compute(qoos::StepQlaOoS, inputs::LabeledValues, withJac::Bool)::MMResult
    # Build labels
    Dls = collect(DkLabel(qoos.material, qoos.shell, k) for k = 1:3)
    Pls = collect(PkLabel(qoos.material, qoos.shell, k) for k = 1:3)
    Tls = collect(TkLabel(qoos.material, qoos.shell, k) for k = 1:3)
    e0l, ml = E0keVLabel(qoos.material), mLabel(qoos.shell)
    Jl, Ml, Zl = JLabel(qoos.material), BigMLabel(qoos.material), ZbarbLabel(qoos.material)
    # Extract values
    D, P, T = map(l -> inputs[l], Dls), map(l -> inputs[l], Pls), map(l -> inputs[l], Tls)
    E0, m, M, J, Zbarb = inputs[e0l], inputs[ml], inputs[Ml], inputs[Jl], inputs[Zl]
    Ea = 0.001 * energy(qoos.shell) # in keV
    U0, V0 = E0 / Ea, E0 / J
    # Calculate results
    Qla = log(U0) / ((U0^m) * Ea^2)
    h = collect((D[k] * (V0 / U0)^P[k]) * (T[k] * U0^T[k] * log(U0) - U0^T[k] + 1.0) / (T[k]^2) for k in eachindex(Dls))
    f = J / (M * Ea)
    OoS = f * sum(h)
    η = 1.75e-3 * Zbarb + 0.37 * (1.0 - exp(-0.015 * Zbarb^1.3))
    JU0 = 1.0 + U0 * (log(U0) - 1.0)
    Wbar = 0.595 + η / 3.7 + η^4.55
    q = (2.0 * Wbar - 1.0) / (1.0 - Wbar)

    Qlal, OoSl = QlaLabel(qoos.material, qoos.shell), OoSLabel(qoos.material, qoos.shell)
    ηl, JU0l = ηLabel(qoos.material), JU0Label(qoos.material)
    Wbarl, ql = WbarLabel(qoos.material), qLabel(qoos.material)
    vals = LabeledValues([Qlal, OoSl, ηl, JU0l, Wbarl, ql], [Qla, OoS, η, JU0, Wbar, q])

    jacob = withJac ? zeros(Float64, length(vals), length(inputs)) : missing
    if withJac
        # Qla (ok!)
        @Assert indexin(Qlal, vals) == 1
        jacob[1, indexin(e0l, inputs)] = (1.0 - m * log(U0)) / (U0^(m + 1) * Ea^3)
        jacob[1, indexin(ml, inputs)] = -Qla * log(U0)
        # 1/S (ok!)
        @Assert indexin(OoSl, vals) == 2
        for k in eachindex(h)
            jacob[2, indexin(Dls[k], inputs)] = f * (h[k] / D[k])
            jacob[2, indexin(Pls[k], inputs)] = f * log(Ea / J) * h[k]
            jacob[2, indexin(Tls[k], inputs)] =
                f *
                (1 / T[k]) * #
                (D[k] * (U0^T[k]) * (Ea / J)^P[k] * log(U0)^2 - 2.0 * h[k])
        end
        jacob[2, indexin(Ml, inputs)] = -OoS / M
        jacob[2, indexin(e0l, inputs)] =
            (f / E0) * log(U0) * sum(D[k] * ((Ea / J)^P[k]) * (U0^T[k]) for k in eachindex(h))
        jacob[2, indexin(Jl, inputs)] = (-1.0 / (M * Ea)) * sum(P[k] * h[k] for k in eachindex(h))
        # η (ok!)
        @Assert indexin(ηl, vals) == 3
        δηδZbarb = 1.75e-3 + 7.215e-3 * exp(-0.015 * Zbarb^1.3) * (Zbarb^0.3)
        jacob[3, indexin(Zl, inputs)] = δηδZbarb
        # J(U0) (ok!)
        @Assert indexin(JU0l, vals) == 4
        jacob[4, indexin(e0l, inputs)] = log(U0) / Ea
        # Wbar (ok!)
        @Assert indexin(Wbarl, vals) == 5
        δWbarδZbarb = (1.0 / 3.7 + 4.55 * η^3.55) * δηδZbarb
        jacob[5, indexin(Zl, inputs)] = δWbarδZbarb
        # q (ok!)
        @Assert indexin(ql, vals) == 6
        jacob[6, indexin(Zl, inputs)] = (1.0 / ((Wbar - 1.0)^2)) * δWbarδZbarb
    end
    return (vals, jacob)
end

struct StepRPhi0 <: MeasurementModel
    material::String
    shell::AtomicSubShell
end

struct RLabel <: Label
    material::String
    shell::AtomicSubShell
end
Base.show(io::IO, rl::RLabel) = print(io, "R[$(rl.material),$(rl.shell)]")

struct ϕ0Label <: Label
    material::String
    shell::AtomicSubShell
end
Base.show(io::IO, l::ϕ0Label) = print(io, "ϕ₀[$(l.material),$(l.shell)]")

function NeXLUncertainties.compute(rp::StepRPhi0, inputs::LabeledValues, withJac::Bool)::MMResult
    # inputs
    E0l, ql, JU0l = E0keVLabel(rp.material), qLabel(rp.material), JU0Label(rp.material)
    ηl, Wbarl = ηLabel(rp.material), WbarLabel(rp.material)
    e0, q, JU0 = inputs[E0l], inputs[ql], inputs[JU0l]
    η, Wbar = inputs[ηl], inputs[Wbarl]
    Ea = 0.001 * energy(rp.shell)
    U0 = e0 / Ea
    # outputs
    GU0 = (U0 - 1.0 - (1.0 - U0^(-1.0 - q)) / (1 + q)) / ((2.0 + q) * JU0)
    R = 1.0 - η * Wbar * (1.0 - GU0)
    ϕ0 = 1.0 + 3.3 * (1.0 - 1.0 / U0^(2.0 - 2.3 * η)) * η^1.2

    Rl, ϕ0l = RLabel(rp.material, rp.shell), ϕ0Label(rp.material, rp.shell)
    vals = LabeledValues([Rl, ϕ0l], [R, ϕ0])
    jac = withJac ? zeros(Float64, length(vals), length(inputs)) : missing
    if withJac
        # δGUδq = ((U0^(1.0+q)-1.0)-(1.0-q)*log(U0))/(JU0*(1.0+q)^2*(2.0+q)*U0^(1.0+q))
        δGUδE0 = (1.0 - U0^(-2.0 - q)) / (Ea * JU0 * (2.0 + q))
        δGUδq = ((1.0 - U0^(-1 - q)) / ((1 + q)^2) - (U0^(-1 - q) * log(U0)) / (1 + q) - GU0 * JU0) / (JU0 * (2 + q))
        δGUδJU0 = -GU0 / JU0
        @Assert indexin(Rl, vals) == 1
        jac[1, indexin(ηl, inputs)] = Wbar * (GU0 - 1.0) # δRδη (ok)
        jac[1, indexin(E0l, inputs)] = η * Wbar * δGUδE0 # δRδE0
        jac[1, indexin(ql, inputs)] = η * Wbar * δGUδq   # δRδq
        jac[1, indexin(JU0l, inputs)] = η * Wbar * δGUδJU0   # δRδU0
        jac[1, indexin(Wbarl, inputs)] = η * (GU0 - 1.0) # δRδWbar
        @Assert indexin(ϕ0l, vals) == 2
        jac[2, indexin(ηl, inputs)] = 3.96 * η^0.2 * (1.0 - U0^(2.3η - 2.0)) - 7.59 * η * U0^(2.3η - 2.0) * log(U0) # δϕ0δη
        jac[2, indexin(E0l, inputs)] = (7.59 * η^1.2 * (0.869565 - η) * U0^(2.3η - 3.0)) / Ea # δϕ0δE0
    end
    return (vals, jac)
end

struct StepFRBar <: MeasurementModel
    material::String
    shell::AtomicSubShell
end

struct RbarLabel <: Label
    material::String
    shell::AtomicSubShell
end
Base.show(io::IO, l::RbarLabel) = print(io, "Rbar[$(l.material),$(l.shell)]")

struct FLabel <: Label
    material::String
    shell::AtomicSubShell
end
Base.show(io::IO, l::FLabel) = print(io, "F[$(l.material),$(l.shell)]")

function NeXLUncertainties.compute(st::StepFRBar, inputs::LabeledValues, withJac::Bool)::MMResult
    # Input labels
    Zbl, Qlal = ZbarbLabel(st.material), QlaLabel(st.material, st.shell)
    OoSl, Rl = OoSLabel(st.material, st.shell), RLabel(st.material, st.shell)
    ϕ0l, e0l = ϕ0Label(st.material, st.shell), E0keVLabel(st.material)
    # Input values
    Zb, Qla, OoS, R, ϕ0 = inputs[Zbl], inputs[Qlal], inputs[OoSl], inputs[Rl], inputs[ϕ0l]
    E0, Ea = inputs[e0l], 0.001 * energy(st.shell)
    U0 = E0 / Ea
    # Computed quantities
    X, Y, u42 = 1.0 + 1.3 * log(Zb), 0.2 + Zb / 200.0, U0^0.42
    FoRbar = 1.0 + (X * log(1.0 + Y * (1.0 - 1.0 / u42))) / log(1.0 + Y)
    F = R * OoS / Qla
    if FoRbar >= ϕ0
        Rbar = F / FoRbar
    else
        Rbar = R / ϕ0
    end
    # Output quantities
    fl, rbarl = FLabel(st.material, st.shell), RbarLabel(st.material, st.shell)
    vals = LabeledValues([fl, rbarl], [F, Rbar])
    jac = withJac ? zeros(Float64, length(vals), length(inputs)) : missing
    if withJac
        δXδZb, δYδZb, = 1.3 / Zb, 1.0 / 200.0
        δFoRbδX = log(1.0 + Y * (1.0 - 1.0 / u42)) / log(1.0 + Y)
        δFoRbδY =
            X * (((u42 - 1.0) * log(1.0 + Y)) / (u42 + Y * (u42 - 1.0)) - log(1.0 + Y - Y / u42) / (1.0 + Y)) /
            (log(1 + Y)^2)

        δFoRbδZb = δFoRbδX * δXδZb + δFoRbδY * δYδZb
        δFoRbδU0 = (0.42 * X * Y) / (U0 * u42 * (1.0 + Y * (1.0 - 1.0 / u42)) * log(1.0 + Y))

        δFδR = OoS / Qla
        δFδOoS = R / Qla
        δFδQla = -R * OoS / (Qla^2)

        if FoRbar >= ϕ0
            δRbarδF = 1.0 / FoRbar
            δRbarδFoRbar = -F / FoRbar^2
            δRbarδϕ0 = 0.0
        else
            δRbarδF = 1.0 / ϕ0
            δRbarδϕ0 = -F / ϕ0^2
            δRbarδFoRbar = 0.0
        end
        @Assert indexin(fl, vals) == 1
        jac[1, indexin(Qlal, inputs)] = δFδQla
        jac[1, indexin(OoSl, inputs)] = δFδOoS
        jac[1, indexin(Rl, inputs)] = δFδR
        jac[1, indexin(Zbl, inputs)] = δRbarδFoRbar * δFoRbδZb
        jac[1, indexin(ϕ0l, inputs)] = 0.0
        jac[1, indexin(e0l, inputs)] = 0.0
        @Assert indexin(rbarl, vals) == 2
        jac[2, indexin(Qlal, inputs)] = δRbarδF * δFδQla
        jac[2, indexin(OoSl, inputs)] = δRbarδF * δFδOoS
        jac[2, indexin(Rl, inputs)] = δRbarδF * δFδR
        jac[2, indexin(Zbl, inputs)] = δRbarδFoRbar * δFoRbδZb
        jac[2, indexin(ϕ0l, inputs)] = δRbarδϕ0
        jac[2, indexin(e0l, inputs)] = δRbarδFoRbar * δFoRbδU0 / Ea
    end
    return (vals, jac)
end

struct StepPb <: MeasurementModel
    material::String
    shell::AtomicSubShell
end

struct PLabel <: Label
    material::String
    shell::AtomicSubShell
end
Base.show(io::IO, l::PLabel) = print(io, "P[$(l.material),$(l.shell)]")


struct bLabel <: Label
    material::String
    shell::AtomicSubShell
end
Base.show(io::IO, l::bLabel) = print(io, "b[$(l.material),$(l.shell)]")


function NeXLUncertainties.compute(st::StepPb, inputs::LabeledValues, withJac::Bool)::MMResult
    # Extract input variables
    args = (st.material, st.shell)
    Rbarl, Fl, ϕ0l = RbarLabel(args...), FLabel(args...), ϕ0Label(args...)
    Zbarl, e0l = ZbarbLabel(st.material), E0keVLabel(st.material)
    Rbar, F, Zbar, ϕ0, e0 = inputs[Rbarl], inputs[Fl], inputs[Zbarl], inputs[ϕ0l], inputs[e0l]
    Ea = 0.001 * energy(st.shell)
    U0 = e0 / Ea
    # Compute outputs
    g = 0.22 * log(4Zbar) * (1.0 - 2.0 * exp(-Zbar * (U0 - 1.0) / 15.0))
    h = 1.0 - 10.0 * (1.0 - 1.0 / (1.0 + U0 / 10.0)) / (Zbar^2)
    b = sqrt(2.0) * (1.0 + sqrt(1.0 - Rbar * ϕ0 / F)) / Rbar
    gh4max = 0.9 * b * (Rbar^2) * (b - 2.0ϕ0 / F)
    if g * (h^4) < gh4max
        P = g * (h^4) * F / (Rbar^2)
    else
        P = gh4max * F / (Rbar^2)
    end
    vals = LabeledValues([bLabel(args...), PLabel(args...)], [b, P])
    jac = withJac ? zeros(Float64, length(vals), length(inputs)) : missing
    if withJac
        δgδZb = g / (Zbar * log(4.0Zbar)) + 0.0293 * exp(-Zbar * (U0 - 1.0) / 15.0) * (U0 - 1.0) * log(4.0Zbar)
        δgδU0 = 0.0293 * exp(-Zbar * (U0 - 1.0) / 15.0) * Zbar * log(4.0Zbar)
        δhδZb = 20.0 * (1.0 - 1.0 / (1.0 + U0 / 10.0)) / (Zbar^3)
        δhδU0 = -1.0 / (((1.0 + U0 / 10.0) * Zbar)^2)
        δbδRbar = -b / Rbar - ϕ0 / (F * Rbar * sqrt(2.0 * (1.0 - ϕ0 * Rbar / F)))
        δbδϕ0 = -1.0 / (F * sqrt(2.0 * (1.0 - ϕ0 * Rbar / F)))
        δbδF = ϕ0 / (F^2 * sqrt(2.0 * (1.0 - ϕ0 * Rbar / F)))

        δPδg, δPδh, δPδF, δPδRbar = P / g, 4.0P / h, P / F, -2.0 * P / Rbar

        @Assert indexin(bLabel(args...), vals) == 1
        jac[1, indexin(Rbarl, inputs)] = δbδRbar
        jac[1, indexin(Fl, inputs)] = δbδF
        jac[1, indexin(ϕ0l, inputs)] = δbδϕ0
        # jac[1, indexin(Zbarl, inputs)] = 0.0
        # jac[1, indexin(e0l, inputs)] = 0.0
        @Assert indexin(PLabel(args...), vals) == 2
        if g * (h^4) < gh4max
            jac[2, indexin(Rbarl, inputs)] = δPδRbar
            jac[2, indexin(Fl, inputs)] = δPδF
            # jac[2, indexin(ϕ0l, inputs)] = 0.0
            jac[2, indexin(Zbarl, inputs)] = δPδg * δgδZb + δPδh * δhδZb
            jac[2, indexin(e0l, inputs)] = (δPδg * δgδU0 + δPδh * δhδU0) / Ea
        else
            δPδϕ0 = ((1.8 * Rbar) / F) * ((F * b - ϕ0) * δbδϕ0 - b)
            #jac[2, indexin(Rbarl, inputs)] = 0.0
            jac[2, indexin(Fl, inputs)] = 0.9 * b^2
            jac[2, indexin(ϕ0l, inputs)] = δPδϕ0
            #jac[2, indexin(Zbarl, inputs)] = 0.0
            #jac[2, indexin(e0l, inputs)] = 0.0
        end
    end
    return (vals, jac)
end

struct Stepaϵ <: MeasurementModel
    material::String
    shell::AtomicSubShell
end

struct aLabel <: Label
    material::String
    shell::AtomicSubShell
end
Base.show(io::IO, l::aLabel) = print(io, "a[$(l.material),$(l.shell)]")

struct ϵLabel <: Label
    material::String
    shell::AtomicSubShell
end
Base.show(io::IO, l::ϵLabel) = print(io, "ϵ[$(l.material),$(l.shell)]")

function NeXLUncertainties.compute(st::Stepaϵ, inputs::LabeledValues, withJac::Bool)::MMResult
    # Extract input variables
    args = (st.material, st.shell)
    Rbarl, Fl, ϕ0l = RbarLabel(args...), FLabel(args...), ϕ0Label(args...)
    Pl, bl = PLabel(args...), bLabel(args...)
    Rbar, F, ϕ0, P, b = inputs[Rbarl], inputs[Fl], inputs[ϕ0l], inputs[Pl], inputs[bl]
    # Compute the values
    den, tiny = b * F * (2.0 - b * Rbar) - ϕ0, 1.0e-6
    a = (P + b * (2.0 * ϕ0 - b * F)) / den
    ϵ = (a - b) / b
    if abs(ϵ) < tiny
        ϵ = tiny
        a = b * (1.0 + ϵ)
    end
    vals = LabeledValues([aLabel(args...), ϵLabel(args...)], [a, ϵ])
    jac = withJac ? zeros(Float64, length(vals), length(inputs)) : missing
    if withJac
        δaδP = 1.0 / den # Ok
        δaδϕ0 = (b^2 * F * (3.0 - 2.0 * b * Rbar) + P) / (den^2) # Ok
        δaδRbar = (b^2 * F * (P + 2.0 * b * ϕ0 - b^2 * F)) / (den^2) # Ok
        δaδF = b * (P * (b * Rbar - 2.0) + b * ϕ0 * (2.0 * b * Rbar - 3.0)) / (den^2) # Ok
        δaδb = -2.0 * (F * (P - b * (ϕ0 + P * Rbar) + b^2 * (F - ϕ0 * Rbar)) + ϕ0^2) / (den^2) # Ok

        δϵδa = 1.0 / b # Ok
        δϵδb = -a / (b^2) # Ok

        @Assert indexin(aLabel(args...), vals) == 1
        jac[1, indexin(Rbarl, inputs)] = δaδRbar
        jac[1, indexin(Fl, inputs)] = δaδF
        jac[1, indexin(ϕ0l, inputs)] = δaδϕ0
        jac[1, indexin(Pl, inputs)] = δaδP
        jac[1, indexin(bl, inputs)] = δaδb
        @Assert indexin(ϵLabel(args...), vals) == 2
        jac[2, indexin(Rbarl, inputs)] = δϵδa * δaδRbar
        jac[2, indexin(Fl, inputs)] = δϵδa * δaδF
        jac[2, indexin(ϕ0l, inputs)] = δϵδa * δaδϕ0
        jac[2, indexin(Pl, inputs)] = δϵδa * δaδP
        jac[2, indexin(bl, inputs)] = (b * δaδb - a) / (b^2)
    end
    return (vals, jac)
end

struct StepAB <: MeasurementModel
    material::String
    shell::AtomicSubShell
end

struct ALabel <: Label
    material::String
    shell::AtomicSubShell
end
Base.show(io::IO, l::ALabel) = print(io, "A[$(l.material),$(l.shell)]")

struct BLabel <: Label
    material::String
    shell::AtomicSubShell
end
Base.show(io::IO, l::BLabel) = print(io, "B[$(l.material),$(l.shell)]")

function NeXLUncertainties.compute(st::StepAB, inputs::LabeledValues, withJac::Bool)::MMResult
    # Build input variable labels
    args = (st.material, st.shell)
    ϵl, bl, Fl, Pl, ϕ0l = ϵLabel(args...), bLabel(args...), FLabel(args...), PLabel(args...), ϕ0Label(args...)
    # Extract input variables
    ϵ, b, F, P, ϕ0 = inputs[ϵl], inputs[bl], inputs[Fl], inputs[Pl], inputs[ϕ0l]
    # Compute the values
    B = (b^2 * F * (1.0 + ϵ) - P - ϕ0 * b * (2.0 + ϵ)) / ϵ
    A = (B / b + ϕ0 - b * F) * ((1.0 + ϵ) / ϵ)
    Al, Bl = ALabel(args...), BLabel(args...)
    vals = LabeledValues([Al, Bl], [A, B])
    jac = withJac ? zeros(Float64, length(vals), length(inputs)) : missing
    if withJac
        δBδϵ = (P + 2.0 * b * ϕ0 - b^2 * F) / (ϵ^2)
        δBδb = (2.0 * b * F * (1.0 + ϵ) - (2.0 + ϵ) * ϕ0) / ϵ
        δBδF = b^2 * (1.0 + ϵ) / ϵ
        δBδP = -1.0 / ϵ
        δBδϕ0 = -b * (2.0 + ϵ) / ϵ

        δAδϵ = -A / (ϵ * (1.0 + ϵ)) + (δBδϵ / b) * ((1.0 + ϵ) / ϵ)
        δAδb = (δBδb / b - (B / (b^2) + F)) * ((1.0 + ϵ) / ϵ)
        δAδF = ((δBδF - b^2) / b) * ((1.0 + ϵ) / ϵ)
        δAδϕ0 = ((b + δBδϕ0) / b) * ((1.0 + ϵ) / ϵ)
        δAδP = (δBδP / b) * ((1.0 + ϵ) / ϵ)

        Ai, Bi = indexin(Al, vals), indexin(Bl, vals)
        jac[Ai, indexin(ϵl, inputs)] = δAδϵ
        jac[Ai, indexin(bl, inputs)] = δAδb
        jac[Ai, indexin(Fl, inputs)] = δAδF
        jac[Ai, indexin(Pl, inputs)] = δAδP
        jac[Ai, indexin(ϕ0l, inputs)] = δAδϕ0

        jac[Bi, indexin(ϵl, inputs)] = δBδϵ
        jac[Bi, indexin(bl, inputs)] = δBδb
        jac[Bi, indexin(Fl, inputs)] = δBδF
        jac[Bi, indexin(Pl, inputs)] = δBδP
        jac[Bi, indexin(ϕ0l, inputs)] = δBδϕ0
    end
    return (vals, jac)
end

struct StepχFr <: MeasurementModel
    material::String
    xray::CharXRay
end

struct χLabel <: Label
    material::String
    xray::CharXRay
end
Base.show(io::IO, l::χLabel) = print(io, "χ[$(l.material),$(l.xray)]")

struct θLabel <: Label
    material::String
end
Base.show(io::IO, l::θLabel) = print(io, "θ[$(l.material)]")

struct dzLabel <: Label
    material::String
end
Base.show(io::IO, l::dzLabel) = print(io, "Δz[$(l.material)]")

struct FrLabel <: Label
    material::String
    xray::CharXRay
end
Base.show(io::IO, l::FrLabel) = print(io, "Fᵣ[$(l.material),$(l.xray)]")

function NeXLUncertainties.compute(st::StepχFr, inputs::LabeledValues, withJac::Bool)::MMResult
    # Build input variable labels
    args = (st.material, inner(st.xray))
    μoρl, θl, bl, Bl = μoρLabel(st.material, st.xray), θLabel(st.material), bLabel(args...), BLabel(args...)
    Al, ϕ0l, ϵl, dzl = ALabel(args...), ϕ0Label(args...), ϵLabel(args...), dzLabel(st.material)
    # Extract input variables
    μoρ, θ, b, B = inputs[μoρl], inputs[θl], inputs[bl], inputs[Bl]
    A, ϕ0, ϵ = inputs[Al], inputs[ϕ0l], inputs[ϵl]
    # Compute the values
    χ = μoρ * csc(θ)
    b1ϵχ, bχ = (b * (1.0 + ϵ) + χ), b + χ
    Fr = (ϕ0 + B / bχ - A * b * ϵ / b1ϵχ) / bχ
    χl, Frl = χLabel(st.material, st.xray), FrLabel(st.material, st.xray)
    vals = LabeledValues([χl, Frl], [χ, Fr])
    jac = withJac ? zeros(Float64, length(vals), length(inputs)) : missing
    if withJac
        δχδμoρl = χ / μoρ
        δχδθ = -χ * cot(θ)
        δFrδχ = (-2.0 * B + (bχ * A * b * ϵ * (2.0 * χ + b * (2.0 + ϵ)) - b1ϵχ^2 * ϕ0) / (b1ϵχ^2)) / (bχ^3)
        χi, Fri = indexin(χl, vals), indexin(Frl, vals)
        jac[χi, indexin(μoρl, inputs)] = δχδμoρl
        jac[χi, indexin(θl, inputs)] = δχδθ
        jac[Fri, indexin(μoρl, inputs)] = δFrδχ * δχδμoρl
        jac[Fri, indexin(θl, inputs)] = δFrδχ * δχδθ
        jac[Fri, indexin(dzl, inputs)] = -χ * Fr
        jac[Fri, indexin(Bl, inputs)] = 1.0 / (bχ^2)
        jac[Fri, indexin(bl, inputs)] = (-B / (bχ^2) + (-A * χ * ϵ) / (b1ϵχ^2) - Fr) / bχ
        jac[Fri, indexin(Al, inputs)] = (-b * ϵ) / (bχ * b1ϵχ)
        jac[Fri, indexin(ϕ0l, inputs)] = 1.0 / bχ
        jac[Fri, indexin(ϵl, inputs)] = (-A * b) / (b1ϵχ^2)
    end
    return (vals, jac)
end

struct StepFrc <: MeasurementModel
    material::String
    coating::String
    xray::CharXRay
end

struct FrcLabel <: Label
    material::String
    coating::String
    xray::CharXRay
end
Base.show(io::IO, l::FrcLabel) = print(io, "Fᵣᶜ[$(l.material),$(l.xray),$(l.coating)]")

struct tcLabel <: Label
    coating::String
end
Base.show(io::IO, l::tcLabel) = print(io, "tᶜ[$(l.coating)]")

function NeXLUncertainties.compute(st::StepFrc, inputs::LabeledValues, withJac::Bool)::MMResult
    # Build input variable labels
    θl, tcl = θLabel(st.material), tcLabel(st.coating)
    Frl, μoρcl = FrLabel(st.material, st.xray), μoρLabel(st.coating, st.xray)
    # Extract input variables
    θ, tc, Fr, μoρc = inputs[θl], inputs[tcl], inputs[Frl], inputs[μoρcl]
    # Compute the values
    Frc = exp(-μoρc * tc * csc(θ)) * Fr
    Frcl = FrcLabel(st.material, st.coating, st.xray)
    vals = LabeledValues([Frcl], [Frc])
    jac = withJac ? zeros(Float64, length(vals), length(inputs)) : missing
    if withJac
        Frci = indexin(Frcl, vals)
        @Assert Frci == 1
        jac[Frci, indexin(tcl, inputs)] = -Frc * μoρc * csc(θ)
        jac[Frci, indexin(θl, inputs)] = Frc * μoρc * tc * csc(θ) * cot(θ)
        jac[Frci, indexin(μoρcl, inputs)] = -Frc * tc * csc(θ)
        jac[Frci, indexin(Frl, inputs)] = Frc / Fr
    end
    return (vals, jac)
end

struct StepZA <: MeasurementModel
    unknown::String
    standard::String
    xray::CharXRay
    coatingU::String
    coatingS::String
end

struct ZLabel <: Label
    unknown::String
    standard::String
    shell::AtomicSubShell
end
Base.show(io::IO, l::ZLabel) = print(io, "Z[$(l.unknown),$(l.standard)]")

struct AbsLabel <: Label
    unknown::String
    standard::String
    xray::CharXRay
    coatingU::String
    coatingS::String
end
Base.show(io::IO, l::AbsLabel) = print(io, "A[$(l.unknown),$(l.standard),$(l.xray)]")

struct ZALabel <: Label
    unknown::String
    standard::String
    xray::CharXRay
    coatingU::String
    coatingS::String
end
Base.show(io::IO, l::ZALabel) = print(io, "ZA[$(l.unknown),$(l.standard),$(l.xray)]")

function NeXLUncertainties.compute(st::StepZA, inputs::LabeledValues, withJac::Bool)::MMResult
    # Build input variable labels
    shell = inner(st.xray)
    Ful, Frcul = FLabel(st.unknown, shell), FrcLabel(st.unknown, st.coatingU, st.xray)
    Fsl, Frcsl = FLabel(st.standard, shell), FrcLabel(st.standard, st.coatingS, st.xray)
    # Extract input variables
    Fu, Frcu, Fs, Frcs = inputs[Ful], inputs[Frcul], inputs[Fsl], inputs[Frcsl]
    # Compute the values
    Z, A, ZA = Fu / Fs, (Frcu / Fu) / (Frcs / Fs), Frcu / Frcs
    Zl = ZLabel(st.unknown, st.standard, shell)
    Al = AbsLabel(st.unknown, st.standard, st.xray, st.coatingU, st.coatingS)
    ZAl = ZALabel(st.unknown, st.standard, st.xray, st.coatingU, st.coatingS)
    vals = LabeledValues([Zl, Al, ZAl], [Z, A, ZA])
    jac = withJac ? zeros(Float64, length(vals), length(inputs)) : missing
    if withJac
        Zi = indexin(Zl, vals)
        jac[Zi, indexin(Ful, inputs)] = 1.0 / Fs
        jac[Zi, indexin(Fsl, inputs)] = -Z / Fs
        Ai = indexin(Al, vals)
        jac[Ai, indexin(Ful, inputs)] = -A / Fu
        jac[Ai, indexin(Fsl, inputs)] = A / Fs
        jac[Ai, indexin(Frcul, inputs)] = A / Frcu
        jac[Ai, indexin(Frcsl, inputs)] = -A / Frcs
        ZAi = indexin(ZAl, vals)
        jac[ZAi, indexin(Frcul, inputs)] = 1.0 / Frcs
        jac[ZAi, indexin(Frcsl, inputs)] = -ZA / Frcs
    end
    return (vals, jac)
end

mjz(sample, elms, all) = StepMJZbarb(sample, elms) | MaintainInputs(mLabel(sample))
dpt(sample, all) = StepDPT(sample, inner(cxr)) | allinp
qla(sample, all) = StepQlaOoS(sample, inner(cxr)) | MaintainInputs([E0keVLabel(sample), ZbarbLabel(sample)])
rp(sample, shell, all) =
    StepRPhi0(sample, shell) |
    MaintainInputs([E0keVLabel(sample), ZbarbLabel(sample), OoSLabel(sample, shell), QlaLabel(sample, shell)], all)
frbar(sample, shell, all) =
    StepFRBar(sample, inner(cxr)) |
    MaintainInputs([ZbarbLabel(sample), ϕ0Label(sample, shell), E0keVLabel(sample)], all)
pb(sample, shell, all) =
    StepPb(sample, shell) |
    MaintainInputs([ϕ0Label(sample, shell), RbarLabel(sample, shell), FLabel(sample, shell)], all)
aϵ(sample, shell, all) =
    Stepaϵ(sample, shell) |
    MaintainInputs([ϕ0Label(sample, shell), FLabel(sample, shell), PLabel(sample, shell), bLabel(sample, shell)], all)
AB(sample, shell, all) =
    NeXLMatrixCorrection.StepAB(unknown, inner(cxr)) |
    MaintainInputs([bLabel(sample, shell), ϕ0Label(sample, shell), ϵLabel(sample, shell), FLabel(sample, shell)], all)

"""
    steps1(sample, elms, shell, all)

steps1 requires as data MassFractionLabel, AtomicWeightLabel, JzLabel, E0Label, mLabel in an UncertainValues
"""
steps1(sample, elms, shell, all = false) =
    AB(sample, shell, all) ∘ aϵ(sample, shell, all) ∘ pb(sample, shell, all) ∘ frbar(sample, shell, all) ∘
    rp(sample, shell, all) ∘ qla(sample, all) ∘ dpt(sample, all) ∘ mjz(sample, elms, all)


χFr(sample, shell, all) =
    NeXLMatrixCorrection.StepχFr(sample, shell) | MaintainInputs([θLabel(sample), FLabel(sample, shell)], all)

"""
    steps2(sample, shell, all)

steps2 requires as data μoρLabel, dzLabel in an UncertainValues
"""
steps2(sample, shell, all = false) = χFr(sample, shell, all)

Frc(sample, xray, coating, all) = StepFrc(sample, coating, xray) | MaintainInputs([FLabel(sample, inner(xray))], all)

"""
    steps3(sample, xray, layer, all)

steps3 requires as data tcLabel, μoρLabel for the coating in an UncertainValues
"""
steps3(sample, xray, layer, all = false) = Frc(sample, xray, layer) | MaintainInputs(FLabel(sample, inner(xray)), all)



function xppu(sample::Material, cxr::CharXRay, coating::Film, e0::UncertainValue, all = false)
    shell = inner(cxr)
    m = m(shell)
    input1 = uvs(
        (MassFractionLabel(sample.name, elm) => uv(sample[elm], 0.01 * sample[elm]) for elm in keys(sample))...,
        (AtomicWeightLabel(sample.name, elm) => uv(a(elm, sample), 0.001 * a(elm, sample)) for elm in keys(sample))...,
        (JzLabel(elm) => Ju(elm) for elm in keys(sample))...,
        NeXLMatrixCorrection.E0Label(unknown) => e0,
        NeXLMatrixCorrection.mLabel(shell) => uv(m, 0.01 * m),
    )
    s1 = steps1(sample.name, [keys(sample)...], shell, all)
    r1 = s1(inputs)

    input2 = uvs(
        μoρLabel(coating.material.name) => mac(coating),
        dzLabel(coating.material.name) => convert(UncertainValue, coating.thickness),
    )
    s2 = steps2(sample.name, shell, coating.material.name, all)
    r2 = s2(cat(r1, input2))

    s3 = steps3(sample.name, shell, coating.material.name, all)
    r3 = s3(r2)
    return r3
end
