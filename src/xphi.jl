using SpecialFunctions: erf, erfc

"""
    ρzx(E0::Float64, Ec::Float64, elm::Element)

X-ray range in g/cm² from Eqn 6 Merlet 1994 or Eqn 7 in Merlet 1995
E0, Ec in keV, Z is the atomic number
"""
function ρzx(E0::Float64, Ec::Float64, elm::Element)::Float64
    @assert E0 > Ec "The beam energy must be larger than the edge energy!"
    @assert E0 > Ec && E0 <= 50.0
    a = 1.845e-6 * evalpoly(E0, (2.6, -0.216, 0.015, -0.000387, 0.00000501)) # Ok
    x =  2.2 - 0.0166 * E0 # Ok
    return (a * (E0^x - Ec^x) * (1.0 + 2.0 / (E0^2))) /
           ((1.078 - 0.015 * z(elm)^0.7)^(1.2 - 0.04 * E0)) # Ok
end

"""
    ρzm(ρzx::Float64, E0::Float64, Ec::Float64, em::Union{Element,Material})

Depth of the peak of the ϕ(ρz) curve in g/cm² (Merlet 1994 eqn 5 or Merlet 1995 eqn 7).
E0, Ec in keV, Z is the atomic number
"""
function ρzm(ρzx::Float64, E0::Float64, Ec::Float64, em::Union{Element,Material})::Float64
    U0, Z = E0 / Ec, z(em)
    # return ρzx * (0.1 + 0.35 * exp(-0.07 * Z)) / (1.0 + 10.0 / (U0^10.0)) # 1994
    return ρzx * (0.1 + 0.35 * exp(-0.07 * Z)) / (1.0 + 1000.0 / (Z * U0^10.0)) # 1995
end


"""
    τ(elm::Element, t::Float64)

Transmission coefficient of Zeller and Ruste 1976 (from Merlet 1994 & 1995)
"""
function τ(elm::Element, ρzm::Float64, ρzx::Float64)::Float64
    Z, t = z(elm), ρzm / ρzx
    τ1, τ2 = (1.0 - t)^(4.65 + 0.0356 * Z), (1.0 - t)^(1.112 + 0.00414 * Z^2)
    return τ1 + 4.65 * (τ2 - τ1) / evalpoly(Z, (3.54, 0.0356, -0.00414)) # Ok
end

mexp(sh::AtomicSubShell) = n(sh) == 1 ? 0.95 : 0.8

"""
    Φm(elm::Element, J::Float64, A::Float64, sh::AtomicSubShell, E0::Float64, ρzm::Float64, ρzx::Float64)

The maximum amplitude of the ϕ(ρz) curve (Located at ρzm).
Z, J, A are the atomic number, ionization potential and atomic weight (mass-fraction averaged except J which is log-averaged)
"""
function Φm(
    elm::Element,
    J::Float64, # Mean ionization potential (keV)
    A::Float64, # Atomic weight
    sh::AtomicSubShell,
    E0::Float64, # in keV
    ρzm::Float64,
    ρzx::Float64,
)::Float64
    Z, m, U0 = z(elm), mexp(sh), E0 / (0.001 * energy(sh)) # Ok
    @assert U0 > 1.0  && U0 < 1000.0
    @assert J > 0.021 && J < 1.0 # keV
    Ud = U0 * (1.0 - 1.03e5 * ρzm * Z / (E0^1.61 * J^0.3 * A)) # Ok
    pZ = -0.25 + 0.0787 * Z^0.3 # Ok
    a = 5.0 * (1.0 - 1.0 / (1.0 + Z)^0.8) * ((E0 / 30.0)^pZ) # Ok
    d1, d2 = 1.0 - m, 7.0 - 4.0 * exp(-0.1 * Z) - m # Ok
    b1 =
        (Ud / U0)^d1 *
        (log(Ud) / log(U0)) *
        (1.0 / d1) *
        (1.0 - (1.0 - 1.0 / (Ud^d1)) / (d1 * log(Ud))) # Ok
    b2 =
        (Ud / U0)^d2 *
        (log(Ud) / log(U0)) *
        (1.0 / d2) *
        (1.0 - (1.0 - 1.0 / (Ud^d2)) / (d2 * log(Ud))) # Ok
    QUd, QU0 = log(Ud) / (Ud^m), log(U0) / (U0^m) # Ok
    return τ(elm, ρzm, ρzx) * (QUd / QU0) +
           a * (0.28 * (1.0 - 0.5 * exp(-0.1 * Z)) * b1 + 0.165 * Z^0.6 * b2) # Ok
end

"""
    Φ0(mat::Material, E0::Float64, Ec::Float64)

The value of the ϕ(ρz) curve at ρz=0.
"""
function Φ0(mat::Material, E0::Float64, Ec::Float64)::Float64
    Z, U0 = z(mat), E0 / Ec
    @assert U0 > 1.0
    d1, d2, d3 = 0.02 * Z, 0.1 * Z, 0.4 * Z # Ok
    b1 = (1.0 / d1) * (1.0 - (1.0 - 1.0 / U0^d1) / (d1 * log(U0))) # Ok
    b2 = (1.0 / d2) * (1.0 - (1.0 - 1.0 / U0^d2) / (d2 * log(U0))) # Ok
    b3 = (1.0 / d3) * (1.0 - (1.0 - 1.0 / U0^d3) / (d3 * log(U0))) # Ok
    pZ = -0.25 + 0.0787 * Z^0.3 # Ok
    a =
        (1.87 * Z) *
        evalpoly(log(Z + 1.0), (0.0, -0.00391, 0.00721, -0.001067)) *
        (E0 / 30.0)^pZ # Ok
    return 1.0 + a * (0.27 * b1 + (1.1 + 5.0 / Z) * (b2 - 1.1 * b3)) # Ok 
end

"""
An implementation of Merlet's XPhi matrix correction algorithm as described in Llovett 2010 and Merlet 1994/1995.
"""
struct XPhi <: MatrixCorrection
    subshell::AtomicSubShell
    material::Material
    E0::Float64 # Beam energy (eV)
    Φm::Float64  # Max amplitude
    ρzm::Float64 # ρz at max amplitude
    α::Float64 # Interior side shape parameter
    β::Float64 # Surface side shape parameter



    function XPhi(
        subshell::AtomicSubShell,
        material::Material,
        E0::Float64, # Beam energy (eV)
        Φm::Float64,  # Max amplitude
        ρzm::Float64, # ρz at max amplitude
        α::Float64, # Interior side shape parameter
        β::Float64 # Surface side shape parameter
    )
        new(subshell, material, E0, Φm, ρzm, α, β)
    end

    """
        XPhi(mat::Material, sh::AtomicSubShell, e0::Float64)

    Construct an object representing the XPhi molde for the specified material and atomic shell at the
    beam energy `e0` in eV.
    """
    function XPhi(mat::Material, sh::AtomicSubShell, e0::Float64)
        E0, Ec = 0.001 * e0, 0.001 * energy(sh)
        mk = keys(mat)
        # This implements the multi-element weighting on page 367 of 1994
        Φ0_v = Φ0(mat, E0, Ec)
        M = sum(nonneg(mat, elm) * z(elm) / a(elm, mat) for elm in mk)
        Φm_v =
            sum(
                nonneg(mat, elm) * z(elm) / a(elm, mat) *
                Φm(elm, 0.001 * J(Berger1982, elm), a(elm, mat), sh, E0,  ρzm(ρzx(E0, Ec, elm), E0, Ec, elm), ρzx(E0, Ec, elm))
                for elm in mk
            ) / M
        ρzx_v = sum(nonneg(mat, elm) * z(elm) / a(elm, mat) * ρzx(E0, Ec, elm) for elm in mk) / M
        ρzm_v = ρzm(ρzx_v, E0, Ec, mat)
        α_v, β_v = 0.46598 * (ρzx_v - ρzm_v), ρzm_v / sqrt(log(Φm_v / Φ0_v))
        @assert isapprox(Φm_v * exp(-(ρzm_v / β_v)^2), Φ0_v, rtol = 1.0e-6)
        @assert isapprox(exp(-((ρzx_v - ρzm_v) / α_v)^2), 0.01, atol = 1.0e-5)
        return new(sh, mat, e0, Φm_v, ρzm_v, α_v, β_v)
    end
end

NeXLCore.range(xp::XPhi) = xp.α * sqrt(-log(0.01)) + xp.ρzm
Base.max(xp::XPhi) = xp.Φm
ϕ0(xp::XPhi) = xp.Φm * exp(-(xp.ρzm / xp.β)^2)

"""
    Φ(xphi::XPhi, ρz::Float64) # W/o absorption
    Φ(xphi::XPhi, ρz::Float64, χ::Float64) # With absorption

The antsatz for ϕ(ρz) curve in the XPhi matrix correction algorithm.
"""
function Φ(xphi::XPhi, ρz::Float64)::Float64
    ρz >= 0.0 ? #
    (
        ρz < xphi.ρzm ? #
        xphi.Φm * exp(-((ρz - xphi.ρzm) / xphi.β)^2) : #
        xphi.Φm * exp(-((ρz - xphi.ρzm) / xphi.α)^2)
    ) : #
    0.0
end
Φ(xphi::XPhi, χ::Float64, ρz::Float64)::Float64 = Φ(xphi, ρz) * exp(-ρz * χ)

"""
    F(xphi::XPhi, 𝒜::Float64, ci::Float64, χ::Float64)

𝒜 - Instrumental, physical and other poorly known parameters
ci - Mass fraction of the i-th element

The generated intensity in the XPhi model.
"""
F(xphi::XPhi, 𝒜::Float64, ci::Float64)::Float64 = #
    0.5 * sqrt(π) * 𝒜 * ci * xphi.Φm * (xphi.α + xphi.β * erf(xphi.ρzm / xphi.β))

"""
    Fχ(xphi::XPhi, 𝒜::Float64, ci::Float64, χ::Float64)


𝒜 - Instrumental, physical and other poorly known parameters
ci - Mass fraction of the i-th element
χ - reduced mass absorption coefficient

The emitted intensity in the XPhi model.
"""
function Fχ(xphi::XPhi, 𝒜::Float64, ci::Float64, χ::Float64)::Float64
    return (0.5 * sqrt(π) * xphi.Φm) * (
        exp(χ * (0.25 * xphi.β^2 * χ - xphi.ρzm)) *
        xphi.β *
        (erf(0.5 * xphi.β * χ) + erf(xphi.ρzm / xphi.β - 0.5 * xphi.β * χ)) +
        exp(χ * (0.25 * xphi.α^2 * χ - xphi.ρzm)) * xphi.α * erfc(0.5 * xphi.α * χ)
    )
end


"""
    Fχ(xphi::XPhi, 𝒜::Float64, ci::Float64, χ::Float64, ρz0::Float64, ρz1::Float64)

𝒜 - Instrumental, physical and other poorly known parameters
ci - Mass fraction of the i-th element
χ - reduced mass absorption coefficient
ρz0, ρz1 - The mass depth range to integrate over

The emitted intensity in the XPhi model for the depth between ρz0 and ρz1
"""
function Fχp(
    xphi::XPhi,
    𝒜::Float64,
    ci::Float64,
    χ::Float64,
    ρz0::Float64,
    ρz1::Float64,
)::Float64
    @assert ρz0 >= 0.0
    @assert ρz0 < ρz1
    if ρz0 < xphi.ρzm && ρz1 > xphi.ρzm
        return Fχp(xphi, 𝒜, ci, χ, ρz0, xphi.ρzm) + Fχp(xphi, 𝒜, ci, χ, xphi.ρzm, ρz1)
    else
        ξ = (ρz0 < xphi.ρzm ? xphi.β : xphi.α)
        return (0.5 * sqrt(π) * xphi.Φm * ξ) *
               exp(χ * (0.25 * ξ^2 * χ - xphi.ρzm)) *
               (
                   erf((ρz1 - xphi.ρzm) / ξ + 0.5 * ξ * χ) -
                   erf((ρz0 - xphi.ρzm) / ξ + 0.5 * ξ * χ)
               )
    end
end


# The remainder of the code is to allow XPhi to fit within the framework with other
# matrix correction algorithms.

matrixcorrection(::Type{XPhi}, mat::Material, ashell::AtomicSubShell, e0::Float64) =
    XPhi(mat, ashell, e0)

function Fχ(xphi::XPhi, xray::CharXRay, θtoa::Float64)::Float64
    @assert inner(xray) == xphi.subshell
    return Fχ(xphi, 1.0, 1.0, χ(material(xphi), xray, θtoa))
end
F(xphi::XPhi) = F(xphi, 1.0, 1.0)
Fχp(xphi::XPhi, xray::CharXRay, θtoa::Real, t::Real) =
    Fχp(xphi, 1.0, 1.0, χ(material(xphi), xray, θtoa), 0.0, t)
ϕ(xphi::XPhi, ρz) = Φ(xphi, ρz)
