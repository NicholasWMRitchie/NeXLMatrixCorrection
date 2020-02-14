using NeXLUncertainties

asserting() = true

# Removeable asserts
macro Assert(test)
  esc(:(if $(@__MODULE__).asserting()
    @assert($test)
   end))
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

Base.show(io::IO, e0::E0keVLabel) = print(io, "E0[", e0.material, " in keV]")

struct ZbarbLabel <: Label
    material::String
end

Base.show(io::IO, zb::ZbarbLabel) = print(io, "Zbarb[", zb.material, "]")

function Ju(elm::Element, f = 0.01) #C1
    j = NeXLMatrixCorrection.J(elm) # from xpp.jl (in eV)
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
    outputs = [BigMLabel(mjz.material), JLabel(mjz.material), ZbarbLabel(mjz.material), E0keVLabel(mjz.material) ]
    M = mapreduce(i -> c[i] * z[i] / a[i], +, eachindex(z))
    J = exp(sum((c[i] * z[i] / a[i]) * log(j[i]) for i in eachindex(z)) / M) # in keV
    Zb = sum(c[i] * sqrt(z[i]) for i in eachindex(z))^2
    values = [M, J, Zb, keV * e0 ]
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

struct DLabel <: Label
    material::String
    shell::AtomicSubShell
    k::Int
end

Base.show(io::IO, d::DLabel) = print(io, "D[", d.material, ",", d.shell, ",", d.k, "]")

struct PLabel <: Label
    material::String
    shell::AtomicSubShell
    k::Int
end

Base.show(io::IO, p::PLabel) = print(io, "P[", p.material, ",", p.shell, ",", p.k, "]")

struct TLabel <: Label
    material::String
    shell::AtomicSubShell
    k::Int
end

Base.show(io::IO, t::TLabel) = print(io, "T[", t.material, ",", t.shell, ",", t.k, "]")

function NeXLUncertainties.compute(dpt::StepDPT, inputs::LabeledValues, withJac::Bool)::MMResult
    # inputs
    Jl, ml = JLabel(dpt.material), mLabel(dpt.shell)
    J, m = inputs[Jl], inputs[ml]
    # outputs
    Dls = [DLabel(dpt.material, dpt.shell, k) for k = 1:3]
    Pls = [PLabel(dpt.material, dpt.shell, k) for k = 1:3]
    Tls = [TLabel(dpt.material, dpt.shell, k) for k = 1:3]

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

struct JU0Label <: Label
    material::String
end

struct WbarLabel <: Label
    material::String
end

struct qLabel <: Label
    material::String
end

struct OoSLabel <: Label
    material::String
    shell::AtomicSubShell
end

Base.show(io::IO, d::OoSLabel) = print(io, "¹/ₛ[", d.material, ",", d.shell, "]")

function NeXLUncertainties.compute(qoos::StepQlaOoS, inputs::LabeledValues, withJac::Bool)::MMResult
    # Build labels
    Dls = collect(DLabel(qoos.material, qoos.shell, k) for k in 1:3)
    Pls = collect(PLabel(qoos.material, qoos.shell, k) for k in 1:3)
    Tls = collect(TLabel(qoos.material, qoos.shell, k) for k in 1:3)
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
        @Assert indexin(Qlal, vals)==1
        jacob[1, indexin(e0l, inputs)] = (1.0 - m * log(U0)) / (U0^(m + 1) * Ea^3)
        jacob[1, indexin(ml, inputs)] = -Qla * log(U0)
        # 1/S (ok!)
        @Assert indexin(OoSl, vals)==2
        for k in eachindex(h)
            jacob[2, indexin(Dls[k], inputs)] = f * (h[k] / D[k])
            jacob[2, indexin(Pls[k], inputs)] = f * log(Ea / J) * h[k]
            jacob[2, indexin(Tls[k], inputs)] = f * (1 / T[k]) * #
                    (D[k] * (U0^T[k]) * (Ea / J)^P[k] * log(U0)^2 - 2.0 * h[k])
        end
        jacob[2, indexin(Ml, inputs)] = -OoS / M
        jacob[2, indexin(e0l, inputs)] = (f / E0) * log(U0) *
                                         sum(D[k] * ((Ea / J)^P[k]) * (U0^T[k]) for k in eachindex(h))
        jacob[2, indexin(Jl, inputs)] = (-1.0 / (M * Ea)) * sum(P[k] * h[k] for k in eachindex(h))
        # η (ok!)
        @Assert indexin(ηl, vals)==3
        δηδZbarb = 1.75e-3 + 7.215e-3 * exp(-0.015 * Zbarb^1.3) * (Zbarb^0.3)
        jacob[3, indexin(Zl, inputs)] = δηδZbarb
        # J(U0) (ok!)
        @Assert indexin(JU0l, vals)==4
        jacob[4, indexin(e0l, inputs)] = log(U0) / Ea
        # Wbar (ok!)
        @Assert indexin(Wbarl, vals)==5
        δWbarδZbarb = (1.0 / 3.7 + 4.55 * η^3.55) * δηδZbarb
        jacob[5, indexin(Zl, inputs)] = δWbarδZbarb
        # q (ok!)
        @Assert indexin(ql, vals)==6
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

Base.show(io::IO, rl::RLabel) = print(io,"R[$(rl.material),$(rl.shell)]")

struct ϕ0Label <: Label
    material::String
    shell::AtomicSubShell
end

Base.show(io::IO, l::ϕ0Label) = print(io,"ϕ₀[$(l.material),$(l.shell)]")

function NeXLUncertainties.compute(rp::StepRPhi0, inputs::LabeledValues, withJac::Bool)::MMResult
    # inputs
    E0l, ql, JU0l = E0keVLabel(rp.material), qLabel(rp.material), JU0Label(rp.material)
    ηl, Wbarl = ηLabel(rp.material), WbarLabel(rp.material)
    e0, q, JU0 = inputs[E0l], inputs[ql], inputs[JU0l]
    η, Wbar = inputs[ηl], inputs[Wbarl]
    Ea = 0.001 * energy(rp.shell)
    U0 = e0/Ea
    # outputs
    GU0 = (U0-1.0-(1.0-U0^(-1.0-q))/(1+q))/((2.0+q)*JU0)
    R = 1.0 - η*Wbar*(1.0-GU0)
    ϕ0 = 1.0+3.3*(1.0-1.0/U0^(2.0-2.3*η))*η^1.2

    Rl, ϕ0l = RLabel(rp.material, rp.shell), ϕ0Label(rp.material, rp.shell)
    vals = LabeledValues( [ Rl, ϕ0l ], [ R, ϕ0 ])
    jac = withJac ? zeros(Float64, length(vals), length(inputs)) : missing
    if withJac
        # δGUδq = ((U0^(1.0+q)-1.0)-(1.0-q)*log(U0))/(JU0*(1.0+q)^2*(2.0+q)*U0^(1.0+q))
        δGUδE0 = (1.0 - U0^(-2.0 - q))/(Ea*JU0*(2.0 + q))
        δGUδq =  ((1.0 - U0^(-1 - q))/((1 + q)^2) - (U0^(-1 - q)*log(U0))/(1 + q) - GU0*JU0)/(JU0*(2 + q))
        δGUδJU0 = -GU0/JU0
        @Assert indexin(Rl, vals) == 1
        jac[1, indexin(ηl, inputs)] = Wbar*(GU0-1.0) # δRδη (ok)
        jac[1, indexin(E0l, inputs)] = η*Wbar*δGUδE0 # δRδE0
        jac[1, indexin(ql, inputs)] = η*Wbar*δGUδq   # δRδq
        jac[1, indexin(JU0l, inputs)] = η*Wbar*δGUδJU0   # δRδU0
        jac[1, indexin(Wbarl, inputs)] = η*(GU0-1.0) # δRδWbar
        @Assert indexin(ϕ0l, vals) == 2
        jac[2, indexin(ηl, inputs)]= 3.96*η^0.2*(1.0-U0^(2.3η-2.0))-7.59*η*U0^(2.3η-2.0)*log(U0) # δϕ0δη
        jac[2, indexin(E0l, inputs)] = (7.59*η^1.2*(0.869565-η)*U0^(2.3η-3.0))/Ea # δϕ0δE0
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

Base.show(io::IO, l::RbarLabel) = print(io,"Rbar[$(l.material),$(l.shell)]")

struct FLabel <: Label
    material::String
    shell::AtomicSubShell
end

Base.show(io::IO, l::FLabel) = print(io,"F[$(l.material),$(l.shell)]")

function NeXLUncertainties.compute(st::StepFRBar, inputs::LabeledValues, withJac::Bool)::MMResult
    # Input labels
    Zbl, Qlal = ZbarbLabel(st.material), QlaLabel(st.material, st.shell)
    OoSl, Rl = OoSLabel(st.material, st.shell), RLabel(st.material, st.shell)
    ϕ0l, e0keV = ϕ0Label(st.material, st.shell), E0keVLabel(st.material)
    # Input values
    Zb, Qla, OoS, R, ϕ0 = inputs[Zbl], inputs[Qlal], inputs[OoSl], inputs[Rl], inputs[ϕ0l]
    E0, Ea = inputs[e0keV], 0.001*energy(st.shell)
    U0 = E0/Ea
    # Computed quantities
    X, Y = 1.0 + 1.3*log(Zb), 0.2 + Zb/200.0
    FoRbar = 1.0 + (X*log(1.0+Y*(1.0-U0^-0.42)))/log(1.0+Y)
    F = R*OoS/Qla
    if FoRbar >= ϕ0
        Rbar = F / FoRbar
    else
        Rbar = R / ϕ0
    end
    # Output quantities
    fl, rbarl =  FLabel(st.material, st.shell), RbarLabel(st.material, st.shell)
    vals = LabeledValues([ fl, rbarl], [F, Rbar])
    jac = withJac ? zeros(Float64, length(vals), length(inputs)) : missing
    if withJac
        δXδZb, δYδZb, = 1.3/Zb, 1.0/200.0
        δFoRbδX = (log(1.0+Y*(1.0-U0^(-0.42))))/log(1.0+Y)
        u42 = U0^0.42
        @show u42
        l1y = log(1.0+Y*(1.0-u42))
        @show l1y
        δFoRbδY = X*((log(1.0+Y)*(u42-1.0))/(u42+Y*(u42-1.0)) - log(1.0+Y*(1.0-u42))/(1.0+Y))/(log(1+Y)^2)
        δFoRbδZb = δFoRbδX*δXδZb  + δFoRbδY*δYδZb
        δFoRbδU0 = (0.42*X*Y)/(U0^-1.42*(1.0+Y*(1.0-U0^-0.42))*log(1.0+Y))

        δFδR = OoS/Qla
        δFδOoS = R/Qla
        δFδQla = -2.0*R*OoS/(Qla^2)
        if FoRbar >= ϕ0
            δRbarδF = 1.0/FoRbar
            δRbarδFoRbar = -F/FoRbar^2
            δRbarδϕ0 = 0.0
        else
            δRbarδF = 1.0/ϕ0
            δRbarδϕ0 = -F/ϕ0^2
            δRbarδFoRbar = 0.0
        end
        @Assert indexin(fl, inputs)==1
        jac[1, Qlal] = δFδQla
        jac[1, OoSl] = δFδOoS
        jac[1, Rl] = δFδR
        jac[1, Zbl] = δRbarδFoRbar * δFoRbδZb
        jac[1, ϕ0l] = 0.0
        jac[1, e0l] = 0.0
        @Assert indexin(rbarl, inputs)==2
        jac[2, Qlal] = δRbarδF * δFδQla
        jac[2, OoSl] = δRbarδF * δFδOoS
        jac[2, Rl] = δRbarδF * δFδR
        jac[2, Zbl] = δRbarδFoRbar * δFoRbδZb
        jac[2, ϕ0l] = δRbarδϕ0
        jac[2, e0l] = δRbarδFoRbar * δFoRbδU0 / Ea
    end
    return (vals, jac)
end
