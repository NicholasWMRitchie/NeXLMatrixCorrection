using SpecialFunctions: erf, erfc

"""
    ρzx1(e0eV::Float64, eceV::Float64, Z::Float64)

X-ray range from Eqn 6 Merlet 1994
"""
function ρzx1(e0eV::Float64, eceV::Float64, Z::Float64)
    # Zero allocations
    E0, Ec = 0.001*e0eV, 0.001*eceV
    a = 1.845e-6*evalpoly(E0, (2.6,-0.216,0.015,-0.000387,0.00000501))
    x = 2.2-0.0166*E0
    return (a*(E0^x-Ec^x)*(1.0+2.0/(E0^2)))/((1.078 - 0.015*Z^0.7)^(1.2-0.04*Ec))
end

"""
    ρzm(ρzx1::Float64, e0eV::Float64, eceV::Float64, Z::Float64)

Position of the maximum of the ϕ(ρz) curve (Merlet 1994 eqn 5).
"""
ρzm(ρzx1::Float64, e0eV::Float64, eceV::Float64, Z::Float64) = #
    ρzx1*(0.1+0.35*exp(-0.07*Z))/(1.0+10.0/((e0eV/eceV)^10.0))
# Zero allocations

"""
    Φm(Z::Float64, sh::AtomicSubShell, e0eV::Float64, ρzm::Float64, ρzx1::Float64)

The maximum value of the ϕ(ρz) curve.
"""
function Φm(Z::Float64, sh::AtomicSubShell, e0eV::Float64, ρzm::Float64, ρzx1::Float64)
    # Zero allocations
    t, E0 = ρzm/ρzx1, 0.001*e0eV
    τ1, τ2 = (1.0 - t)^(4.65+0.0356*Z), (1.0 - t)^(1.112+0.00414*Z^2)
    τ = τ1 + 4.65*(τ2 - τ1)/evalpoly(Z,(3.54,0.0356,0.00414))
    m, U0 = (n(sh)==1 ? 0.95 : 0.7), e0eV/energy(sh)
    j, A = J(Berger1982, element(sh)), NeXLCore.a(element(sh))
    Ud = U0*(1.0-1.03e5*ρzm*Z/(E0^1.61*j^0.3*A))
    pZ = -0.25 + 0.0787*Z^0.3
    a = 5.0*(1.0 - 1.0/(1.0+Z)^0.8)/((E0/30.0)^pZ)
    d1, d2 = 1.0 - m, 7.0 - 4.0*exp(-0.1*Z) - m
    b1 = (Ud/U0)^d1*(log(Ud)/log(U0))*(1.0/d1)*(1.0 - (1.0 - 1.0/(Ud^d1))/(d1*log(Ud)))
    b2 = (Ud/U0)^d2*(log(Ud)/log(U0))*(1.0/d2)*(1.0 - (1.0 - 1.0/(Ud^d2))/(d2*log(Ud)))
    QUd, QU0 = log(Ud)/(Ud^m), log(U0)/(U0^m)
    return τ*(QUd/QU0)+a*(0.28*(1.0-0.5*exp(-0.1*Z))*b1 + 0.165*Z^0.6*b2)
end

"""
    Φ0(Z::Float64, sh::AtomicSubShell, e0eV::Float64)

The value of the ϕ(ρz) curve at ρz=0.
"""
function Φ0(Z::Float64, sh::AtomicSubShell, e0eV::Float64)::Float64
    # Zero allocations
    U0, E0 = e0eV/energy(sh), 0.001*e0eV
    d1, d2, d3 = 0.02*Z, 0.1*Z, 0.4*Z
    b1 = (1.0/d1)*(1.0 - (1.0 - 1.0/U0^d1)/(d1*log(U0)))
    b2 = (1.0/d2)*(1.0 - (1.0 - 1.0/U0^d2)/(d2*log(U0)))
    b3 = (1.0/d3)*(1.0 - (1.0 - 1.0/U0^d3)/(d3*log(U0)))
    pZ = -0.25 + 0.0787*Z^0.3
    a = 1.87*Z*evalpoly(log(Z+1.0),(0.0, -0.00391, 0.00721, -0.001067))*(E0/30.0)^pZ
    return 1.0 + a*(0.27*b1 + (1.1 + 5.0/Z)*(b2-1.1*b3))
end

"""
An implementation of Merlet's XFilm matrix correction algorithm as described in Llovett 2010 and Merlet 1994.
"""
struct XPhi <: MatrixCorrection
    subshell::AtomicSubShell
    material::Material
    E0::Float64 # Beam energy
    Φm::Float64  # Max amplitude
    ρzm::Float64 # ρz at max amplitude
    α::Float64 # Interior side shape parameter
    β::Float64 # Surface side shape parameter

    function XPhi(mat::Material, sh::AtomicSubShell, e0::Float64)
        Z, ec = z(mat), energy(sh)
        ρzx1_v = ρzx1(e0, ec, Z)
        ρzm_v = ρzm(ρzx1_v, e0, ec, Z)
        Φm_v, Φ0_v = Φm(Z, sh, e0, ρzm_v, ρzx1_v), Φ0(Z, sh, e0)
        α_v, β_v = 0.46598*(ρzx1_v - ρzm_v), ρzm_v/sqrt(log(Φm_v/Φ0_v))
        return new(sh, mat, e0, Φm_v, ρzm_v, α_v, β_v)
    end
end

"""
    Φ(xphi::XPhi, ρz::Float64) # W/o absorption
    Φ(xphi::XPhi, ρz::Float64, χ::Float64) # With absorption

The antsatz for ϕ(ρz) curve in the XPhi matrix correction algorithm.
"""
Φ(xphi::XPhi, ρz::Float64) =
    ρz >= 0.0 ? #
        (ρz < xphi.ρzm ? #
            xphi.Φm*exp(-((ρz-xphi.ρzm)/xphi.β)^2) : #
            xphi.Φm*exp(-((ρz-xphi.ρzm)/xphi.α)^2)) : #
        0.0
Φ(xphi::XPhi, χ::Float64, ρz::Float64) =
    Φ(xphi, ρz)*exp(-ρz*χ)

"""
    F(xphi::XPhi, 𝒜::Float64, ci::Float64, χ::Float64)

𝒜 - Instrumental, physical and other poorly known parameters
ci - Mass fraction of the i-th element

The generated intensity in the XPhi model.
"""
F(xphi::XPhi, 𝒜::Float64, ci::Float64) = #
    0.5*sqrt(π)*𝒜*ci*xphi.Φm*(xphi.α + xphi.β*erf(xphi.ρzm/xphi.β))

"""
    Fχ(xphi::XPhi, 𝒜::Float64, ci::Float64, χ::Float64)


𝒜 - Instrumental, physical and other poorly known parameters
ci - Mass fraction of the i-th element
χ - reduced mass absorption coefficient

The emitted intensity in the XPhi model.
"""
Fχ(xphi::XPhi, 𝒜::Float64, ci::Float64, χ::Float64) = #
    0.5*sqrt(π)*xphi.Φm*(
        exp(χ*(0.25*xphi.β^2*χ-xphi.ρzm))*xphi.β*(erf(0.5*xphi.β*χ)+erf(xphi.ρzm/xphi.β-0.5*xphi.β*χ))+
        exp(χ*(0.25*xphi.α^2*χ-xphi.ρzm))*xphi.α*erfc(0.5*xphi.α*χ)
    )

function Fχp(xphi::XPhi, 𝒜::Float64, ci::Float64, χ::Float64, ρz0::Float64, ρz1::Float64)
    res = 0.0
    if ρz0 < xphi.ρzm
        ρmax = min(ρz1, xphi.ρzm)
        res += exp(χ*(0.25*xphi.β^2*χ - xphi.ρzm))*xphi.β*
            (erf((ρ1-xphi.ρzm)/xphi.β+0.5*xphi.β*χ)-erf((ρ0-xphi.ρzm)/xphi.β+0.5*xphi.β*χ))
    end
    if ρz1 > xphi.ρzm
        ρmin = max(ρz0, xphi.ρzm)
        res += exp(χ*(0.25*xphi.α^2*χ - xphi.ρzm))*xphi.α*
            (erf((ρ1-xphi.ρzm)/xphi.α+0.5*xphi.α*χ)-erf((ρ0-xphi.ρzm)/xphi.α+0.5*xphi.α*χ))
    end
    return 0.5*sqrt(π)*xphi.Φm*res
end


matrixcorrection(::Type{XPhi}, mat::Material, ashell::AtomicSubShell, e0::Float64) = XPhi(mat, ashell, e0)

function Fχ(xphi::XPhi, xray::CharXRay, θtoa::Float64)
    @assert inner(xray) == xphi.subshell
    return Fχ(xphi, 1.0, 1.0, χ(material(xphi), xray, θtoa))
end
F(xphi::XPhi) = F(xphi, 1.0, 1.0)
Fχp(xphi::XPhi, xray::CharXRay, θtoa::Real, τ::Real) = #
    Fχp(xphi, 1.0, 1.0, mac(mc.material, xray), 0.0, τ)
ϕ(xphi::XPhi, ρz) = Φ(xphi, ρz)