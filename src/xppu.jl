using NeXLUncertainties

struct JzLabel <: Label
    element::Element
end

struct StepMJZbarb <: MeasurementModel
    material::String
    elements::Vector{Element}
end

struct BigMLabel <: Label
    material::String
end

struct JLabel <: Label
    material::String
end

struct ZbarbLabel <: Label
    material::String
end

function Ju(elm::Element, f=0.01) #C1
    j = NeXLMatrixCorrection.J(elm) # from xpp.jl
    return uv(j,f*j)
end

function NeXLUncertainties.compute(mjz::StepMJZbarb, inputs::LabeledValues, withJac::Bool)::MMResult
    # Build the labels once...
    mfls = map(elm -> MassFractionLabel(mjz.material, elm), mjz.elements)
    awls = map(elm -> AtomicWeightLabel(mjz.material, elm), mjz.elements)
    jzs = map(elm -> JzLabel(elm), mjz.elements)
    # Unpack all the variables
    c, a, = map(mfl -> inputs[mfl], mfls), map(awl -> inputs[awl], awls)
    z, j = map(elm -> convert(Float64, elm.number), mjz.elements),
        map(jz -> inputs[jz], jzs)
    values, jacob = zeros(Float64, 3), withJac ? zeros(Float64, 3, length(inputs)) : missing
    Mi, Ji, Zbi = 1, 2, 3
    outputs = [ BigMLabel(mjz.material), JLabel(mjz.material), ZbarbLabel(mjz.material) ]
    M = (values[Mi] = mapreduce(i -> c[i] * z[i] / a[i], +, eachindex(z)))
    J = (values[Ji] = exp(sum((c[i] * z[i] / a[i]) * log(j[i]) for i in eachindex(z)) / M))
    Zb = (values[Zbi] = sum(c[i] * sqrt(z[i]) for i in eachindex(z))^2)
    if withJac
        for i in eachindex(z)
            mfl, awl, jz = mfls[i], awls[i], jzs[i]
            # dM/d?
            jacob[Mi, indexin(mfl, inputs)] = z[i] / a[i]
            jacob[Mi, indexin(awl, inputs)] = -c[i] * z[i] / (a[i]^2)
            # dJ/d?
            kk = (j[i]*z[i]) / (M*a[i])
            jacob[Ji, indexin(jz, inputs)] = kk * (c[i] / j[i])
            jacob[Ji, indexin(mfl, inputs)] = kk * (log(j[i]) - log(J))
            jacob[Ji, indexin(awl, inputs)] = kk * (c[i] / a[i]) * (log(J) - log(j[i]))
            # dZ/d?
            jacob[Zbi, indexin(mfl, inputs)] = 2.0 * sqrt(z[i] * Zb)
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

struct DLabel <: Label
    material::String
    shell::AtomicSubShell
    k::Int
end

struct PLabel <: Label
    material::String
    shell::AtomicSubShell
    k::Int
end

struct TLabel <: Label
    material::String
    shell::AtomicSubShell
    k::Int
end

function NeXLUncertainties.compute(dpt::StepDPT, inputs::LabeledValues, withJac::Bool)::MMResult
    Dls = ( DLabel(dpt.material,dpt.shell,k) for k in 1:3 )
    Pls = ( PLabel(dpt.material,dpt.shell,k) for k in 1:3 )
    Tls = ( TLabel(dpt.material,dpt.shell,k) for k in 1:3 )
    Jl, ml = JLabel(dpt.material), mLabel(dpt.shell)
    J, m = inputs[Jl], inputs[ml]

    D = ( 6.6e-6, 1.12e-5*(1.35-0.45*J^2), 2.2e-6/J )
    P = ( 0.78, 0.1, 0.25J - 0.5 )
    T = (1.0 - m) .+ P
    lv = LabeledValues([ Dls..., Pls..., Tls...], [ D..., P..., T... ])
    jacob = withJac ? zeros(Float64, length(lv), length(inputs)) : missing
    if withJac
        Jli, mli = indexin(Jl,inputs), indexin(ml, inputs)
        jacob[2, Jli] = -1.008e-5*J # indexin(Dls[2],lv)==2
        jacob[3, Jli] = -2.2e-6/(J^2) # indexin(Dls[3],lv)==3

        jacob[3+3, Jli] = 0.25 # indexin(Pls[3],lv) == 3+3
        jacob[6+3, Jli] = 0.25 # indexin(Tls[3],lv) == 6+3
        for k in 1:3
            jacob[6+k, mli] = -1.0 # indexin(Tls[k],lv)==6+k
        end
    end
    return ( lv, jacob )
end

struct E0Label <: Label
    material::String
end

struct StepQlaOoS <: MeasurementModel
    material::String
    shell::AtomicSubShell
end

struct QlaLabel <: Label
    material::String
    shell::AtomicSubShell
end

struct OoSLabel <: Label
    material::String
    shell::AtomicSubShell
end

function NeXLUncertainties.compute(qoos::StepQlaOoS, inputs::LabeledValues, withJac::Bool)::MMResult
    # Build labels
    Dls = collect( DLabel(qoos.material,qoos.shell,k) for k in 1:3 )
    Pls = collect( PLabel(qoos.material,qoos.shell,k) for k in 1:3 )
    Tls = collect( TLabel(qoos.material,qoos.shell,k) for k in 1:3 )
    e0l, ml = E0Label(qoos.material), mLabel(qoos.shell)
    Jl, Ml = JLabel(qoos.material), BigMLabel(qoos.material)
    # Extract values
    D, P, T = map(l->inputs[l], Dls), map(l->inputs[l], Pls), map(l->inputs[l], Tls)
    E0, m, M, J = inputs[e0l], inputs[ml], inputs[Ml], inputs[Jl]
    Ea = energy(qoos.shell)
    U0, V0 = E0/Ea, E0/J
    # Calculate results
    Qla = log(U0)/((U0^m)*(Ea^2))
    h = collect( (D[k]*(V0/U0)^P[k])*(T[k]*U0^T[k]*log(U0)-U0^T[k]+1.0)/(T[k]^2) for k in 1:3 )
    f = J/(M*Ea)
    OoS = f*sum(h)
    Qlal, OoSl = QlaLabel(qoos.material, qoos.shell), OoSLabel(qoos.material, qoos.shell)
    lv = LabeledValues( [ Qlal, OoSl ], [ Qla, OoS ])

    jacob = withJac ? zeros(Float64, 2, length(inputs)) : missing
    if withJac
        jacob[1, indexin(e0l, inputs)] = (1.0 - m*log(U0))/(U0^(m+1)*Ea^3)
        jacob[1, indexin(ml, inputs)] = -Qla * log(U0)

        for k in 1:3
            jacob[2, indexin(Dls[k], inputs)] = f*(h[k]/D[k])
            jacob[2, indexin(Pls[k], inputs)] = f*log(Ea/J)*h[k]
            jacob[2, indexin(Tls[k], inputs)] = f*((1/T[k])*(D[k]*(E0/Ea)^T[k]*(Ea/J)^P[k])*log(E0/Ea)^2-2.0*h[k])
        end
        jacob[2, indexin(Ml, inputs)] = -OoS/M
        jacob[2, indexin(e0l, inputs)] = (f/E0)*log(E0/Ea)*sum(D[k]*(Ea/J)^P[k]*(E0/Ea)^T[k] for k in 1:3)
        jacob[2, indexin(Jl, inputs)] = (-1.0/(M*Ea))*sum(P[k]*h[k] for k in 1:3)
    end
    return ( lv, jacob )
end
