using PeriodicTable
# Implements Pouchou and Pichoir's XPP model from the description in the book
# Electron Probe Quantification edited by Kurt Heinrich and Dale E. Newbury

"""
    M(mat::Material)
P&P's M from top line in page 35
"""
M(mat::Material) = #C1
    sum(massfraction(elm,mat)*z(elm)/a(elm, mat) for elm in keys(mat))

"""
    J(elm)
Mean ionization potential for the specified element in keV. (PAP Eqn 7)
"""
J(elm::Element) = #C1
    1.0e-3*z(elm)*(10.04+8.25*exp(-z(elm)/11.22))

"""
    J(mat::Material)
Mean ionization potential for the specified material in keV. (PAP Eqn 6)
"""
J(mat::Material) = #C1
    exp(sum(massfraction(elm,mat)*(z(elm)/a(elm,mat))*log(J(elm)) for elm in keys(mat))/M(mat))

# From PAP eqn 8
D(J) = ( 6.6e-6, 1.12e-5*(1.35-0.45*J^2), 2.2e-6/J ) #C1

# From PAP eqn 8
P(J) = ( 0.78, 0.1, -0.5+0.25*J ) #C1

# From PAP near Eqn 11
T(P,m) = collect(1.0 + P[k] - m for k in 1:3) #C1

"""
    dEdρs(args::Dict{Label,AbstractFloat}, mat::MaterialLabel, Ekev::AbstractFloat, elms)
The function P&P uses to describe the deceleration of an electron in a material.
Output units are (keV/cm)/(g/cm^3) = keV cm^2/g. (PAP eqn 5)
"""
function dEdρs(mat::Material, Ekev::AbstractFloat) #C1
    # PAP Eqn 8
    function f(mat::Material, Ekev::AbstractFloat)
        j = J(mat)
        d, p, v = D(j), P(j), Ekev/j
        return sum(i -> d[i]*v^p[i],1:3)
    end
    return -(M(mat)/J(mat))/f(mat,Ekev)
end

"""
    R0(J, D, P, M, Ekev)
Total trajectory (range) of an electron with initial energy Ekev. (in cm/(g/cm^3))
"""
R0(J, D, P, M, Ekev) = #C1
    sum(J^(1.0-P[k])*D[k]*Ekev^(1.0+P[k])/(1.0+P[k]) for k in 1:3)/M


"""
    range(mat::MaterialLabel, e0keV)
Total trajectory (range) of an electron with initial energy e0kev. (in cm/(g/cm^3))
"""
function Base.range(mat::Material, e0keV)
    j = J(mat)
    return R0(j, D(j), P(j), M(mat), e0keV)
end

"""
    ϕ(ρz, A, a, B, b, ϕ0)
Compute the shape of the ϕ(ρz) curve in the XPP model.
"""
ϕ(ρz, A, a, B, b, ϕ0) = #C1
    A*exp(-a*ρz)+(B*ρz+ϕ0-A)*exp(-b*ρz)

"""
    invS(U0, V0, M, D, P, T)
Computes 1/S where S is the stopping power.
"""
invS(U0, V0, M, D, P, T) = #C1
    U0/(V0*M)*sum(D[k]*((V0/U0)^P[k])*(T[k]*U0^T[k]*log(U0)-U0^T[k]+1.0)/(T[k]^2) for k in 1:3)

"""
    S(mat, ashell, E)
Computes S, the stopping power at the specified energy (in keV)
"""
function S(mat::Material, ashell::AtomicShell, Ekev)
    j = J(mat)
    return 1.0 / invS(Ekev/eV2keV(energy(ashell)), Ekev/j, M(mat), D(j), P(j), T(P(j),m(ashell)))
end

"""
    QlA(U,El,m)
Computes the relative ionization cross-section.
"""
QlA(U, El, m) = log(U)/((U^m)*(El^2)) #C1


"""
    m(ashell::AtomicShell)
Returns the ionization cross-section exponent for QlA(U,El,m(ashell))
"""
function m(ashell::AtomicShell) #C1
    if isequal(family(ashell),"K")
        return  0.86+0.12*exp(-(ashell.z/5.0)^2)
    else
        return isequal(family(ashell),"L") ? 0.82 : 0.78
    end
end

# PAP appendix 1
Zbarb(mat::Material) = #C1
    sum(massfraction(elm,mat)*sqrt(z(elm)) for elm in keys(mat))^2

# PAP appendix 1
ηbar(Zbarb) = #C2
    1.75e-3*Zbarb + 0.37*(1.0-exp(-0.015*Zbarb^1.3))

# PAP appendix 1
Wbar(ηbar) = #C2
    0.595 + ηbar/3.7 + ηbar^4.55

q(Wbar) = #C1
    (2.0*Wbar-1.0) / (1.0-Wbar)

# PAP appendix 1
JU0(U0) = #C2
    1.0 + U0*(log(U0)-1.0)

# PAP appendix 1
G(U0, q) = #C2
    (U0 - 1.0 - (1.0 - 1.0/U0^(1.0+q))/(1.0+q))/((2.0+q)*JU0(U0))

# PAP appendix 1
R(ηbar, Wbar, U0) = #C2
    1.0 - ηbar*Wbar*(1.0-G(U0, q(Wbar)))

"""
    R(mat::Material, u0)
Backscatter factor as a function of material and overvoltage.
Compares well to PAP Figure 23.
"""
function R(mat::Material, u0)
    zbarbv = Zbarb(mat)
    ηbarv = ηbar(zbarbv)
    wbar = Wbar(ηbarv)
    return R(ηbarv, wbar, u0)
end

"""
   ϕ0(U0, ηbar)
The value of the ϕ(ρz) curve at ρz=0.0.
"""
ϕ0(U0, ηbar) = #C1
    1.0+3.3*(1.0-1.0/(U0^(2.0-2.3*ηbar)))*ηbar^1.2

X(Zbarb) = #C3
    1.0 + 1.3*log(Zbarb)

Y(Zbarb) = #C3
    0.2 + Zbarb / 200.0

FoverRbar(X, Y, U0) = #C3
    1.0 + X*log(1.0+Y*(1.0-1.0/(U0^0.42))) / log(1.0+Y)

# PAP eqn 13
F(R, invS, QlA) = #C1
    R*invS/QlA

# PAP eqn 28
Rbar(F, FoverRbar, ϕ0) = # C1
    FoverRbar >= ϕ0 ? F / FoverRbar : F / ϕ0

# PAP eqn 28
g(Zbarb, U0) = #C1
    0.22*log(4.0*Zbarb)*(1.0 - 2.0*exp(-Zbarb*(U0-1.0)/15.0))

# PAP eqn 28
h(Zbarb, U0) = #C1
    1.0 - 10.0*(1.0 - 1.0/(1.0+U0/10.0))/(Zbarb^2)

"""
    b(Rbar, ϕ0, F)
XPP ϕ(ρz) model parameter b
""" # PAP Appendix 4
b(Rbar, ϕ0, F) = #C1
    sqrt(2.0)*(1.0+sqrt(1.0-Rbar*ϕ0/F))/Rbar

# PAP eqn 29
P(Zbarb, U0, F, Rbar, ϕ0, b) = #C1
    min(g(Zbarb, U0)*h(Zbarb, U0)^4, 0.9*b*Rbar^2*(b-2.0*ϕ0/F))*F/(Rbar^2)

function ϵ(P, b, ϕ0, F, Rbar) #C1
    a = (P+b*(2.0*ϕ0-b*F))/(b*F*(2.0-b*Rbar)-ϕ0)
    tmp = (a-b)/b
    # print("a = ",a,"(a-b)/b = ",tmp)
    return abs(tmp) > 1.0e-6 ? tmp : sign(tmp)*1.0e-6
end

"""
    a(b, ϵ)
XPP ϕ(ρz) model parameter a
"""
NeXL.a(b, ϵ) = #C1
    b*(1.0+ϵ)

"""
    B(b, F, P, ϕ0, ϵ)
XPP ϕ(ρz) model parameter B
"""
B(b, F, P, ϕ0, ϵ) = #C1
    (b^2*F*(1.0+ϵ) - P - ϕ0*b*(2.0+ϵ))/ϵ

"""
    A(B, b, ϕ0, F, ϵ)
XPP ϕ(ρz) model parameter A
"""
A(B, b, ϕ0, F, ϵ) = #C1
    (B/b + ϕ0 - b*F)*(1+ϵ)/ϵ

function Fχ(χ, A, a, B, b, ϕ0)
    ϵ = (a-b)/b
    (ϕ0 + B / (b+χ) - A*b*ϵ/(b*(1.0+ϵ) + χ))/(b + χ)
end

"""
   XPPCorrection
Represents the essential intermediary values for an XPP matrix correction of
characteristic X-rays from a particular atomic shell in a particular material.
"""
struct XPPCorrection <: MatrixCorrection
    shell::AtomicShell
    material::Material
    E0 # Beam energy (keV)
    ϕ0 # ϕ(ρz) at the surface
    A # Amplitude  factor
    a # Width factor
    B # Amplitude  factor
    b # Width factor
    F # Integral under the ϕ(ρz) curve
    θ # Take-off angle (radiams)
end

Base.show(io::IO, xpp::XPPCorrection) =
    print(io,"XPP[",xpp.shell," in ",name(xpp.material)," at ",xpp.E0," keV, ",rad2deg(xpp.θ),"° TOA]")

"""
    XPP(mat::Material, ashell::AtomicShell, E0)
Construct an XPPCorrection object for the specified material, atomicshell,
beam energy (in keV) and take-off angle (in radians).
"""
function XPP(mat::Material, ashell::AtomicShell, E0, θtoa)
    @assert( (θtoa > 0.0) && (θtoa <= π/2.0), "Take off angle must be in radians in range (0, π/2]")
    Mv, Jv, Zbarbv, eLv = M(mat), J(mat), Zbarb(mat), NeXL.eV2keV(energy(ashell))
    U0v, V0v, ηbarv = E0/eLv, E0/Jv, ηbar(Zbarbv)
    @assert( U0v > 1.0, "The beam energy must be larger than the shell edge energy.")
    @assert( V0v > 1.0, "The beam energy must be larger than the mean energy loss.")
    Dv, Pjv, mv, Wbarv = D(Jv), P(Jv), m(ashell), Wbar(ηbarv)
    invSv = invS(U0v, V0v, Mv, Dv, Pjv, T(Pjv, mv))
    qv, Rv, ϕ0v, QlAv = q(Wbarv), R(ηbarv, Wbarv, U0v), ϕ0(U0v, ηbarv), QlA(U0v, eLv, mv)
    Fv, FoverRbarv = F(Rv, invSv, QlAv), FoverRbar(X(Zbarbv), Y(Zbarbv), U0v)
    Rbarv = Rbar(Fv, FoverRbarv, ϕ0v)
    bv = b(Rbarv, ϕ0v, Fv)
    Pv = P(Zbarbv, U0v, Fv, Rbarv, ϕ0v, bv)
    ϵv = ϵ(Pv, bv, ϕ0v, Fv, Rbarv)
    av, Bv = a(bv, ϵv), B(bv, Fv, Pv, ϕ0v, ϵv)
    Av = A(Bv, bv, ϕ0v, Fv, ϵv)
    return XPPCorrection(ashell, mat, E0, ϕ0v, Av, av, Bv, bv, Fv, θtoa)
end

"""
    xpp(mat::Material, ashell::AtomicShell, E0, θdeg)
Construct an XPPCorrection object for the specified material, atomicshell,
beam energy (in keV) and take-off angle (in degrees).
"""
xpp(mat::Material, ashell::AtomicShell, E0, θtoa)::XPPCorrection =
    XPP(mat, ashell, E0, deg2rad(θtoa))

Fχ(xpp::XPPCorrection, xray::CharXRay) =
    Fχ(χ(material(xpp), xray, takeOffAngle(xpp)), xpp.A, xpp.a, xpp.B, xpp.b, xpp.ϕ0)

F(xpp::XPPCorrection) = xpp.F
atomicShell(mc::XPPCorrection) = mc.shell
takeOffAngle(mc::XPPCorrection) = mc.θ
material(mc::XPPCorrection) = mc.material
beamEnergy(mc::XPPCorrection) = mc.E0

"""
    ϕ(ρz, xpp::XPP)
Computes the ϕ(ρz) curve according to the XPP algorithm.
"""
ϕ(ρz, xpp::XPPCorrection) =
    ϕ(ρz, xpp.A, xpp.a, xpp.B, xpp.b, xpp.ϕ0)

"""
    ϕabs(ρz, xpp::XPP, xray::CharXRay, θtoa)
Computes the absorbed ϕ(ρz) curve according to the XPP algorithm.
"""
ϕabs(ρz, xpp::XPPCorrection, xray::CharXRay) =
    ϕ(ρz, xpp)*exp(-χ(material(xpp), xray, takeOffAngle(xpp))*ρz)

""""
    xppZAF(mat::Material, ashell::AtomicShell, e0, θtoa, coating=NullCoating())
Constructs an ZAFCorrection object using XPPCorrection for the specified parameters.
"""
xppZAF(mat::Material, ashell::AtomicShell, e0, θtoa, coating=NullCoating()) =
    ZAFCorrection(xpp(mat,ashell,e0,θtoa), NullFluorescence(mat, ashell, e0, θtoa), coating)

"""
    buildMultiXPP(mat::Material, cxrs, e0, θtoa, coating=NeXL.NullCoating())
Constructs a MultiZAF around the XPPCorrection algorithm.
"""
function buildMultiXPP(mat::Material, cxrs, e0, θtoa, coating=NeXL.NullCoating())
    zaf(sh) = xppZAF(mat, sh, e0, θtoa, coating)
    zafs = Dict( ( sh, zaf(sh) ) for sh in union(inner.(cxrs)) )
    return MultiZAF(cxrs, zafs)
end
