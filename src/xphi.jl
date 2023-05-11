using SpecialFunctions: erf, erfc

"""
An implementation of Merlet's XPhi matrix correction algorithm as described in Llovett 2010 and Merlet 1994/1995.
"""
struct XPhi <: MatrixCorrection
    subshell::AtomicSubShell
    material::Material{Float64,Float64}
    E0::Float64 # Beam energy (eV)
    Φm::Float64  # Max amplitude
    ρzm::Float64 # ρz at max amplitude
    α::Float64 # Interior side shape parameter
    β::Float64 # Surface side shape parameter



    function XPhi(
        subshell::AtomicSubShell,
        material::Material,
        E0::AbstractFloat, # Beam energy (eV)
        Φm::AbstractFloat,  # Max amplitude
        ρzm::AbstractFloat, # ρz at max amplitude
        α::AbstractFloat, # Interior side shape parameter
        β::AbstractFloat # Surface side shape parameter
    )
        new(subshell, material, E0, Φm, ρzm, α, β)
    end

    """
        XPhi(mat::Material{Float64,Float64}, sh::AtomicSubShell, e0::AbstractFloat)

    Construct an object representing the XPhi molde for the specified material and atomic shell at the
    beam energy `e0` in eV.
    """
    function XPhi(matp::Material, sh::AtomicSubShell, e0::AbstractFloat)
        mat = convert(Material{Float64,Float64}, matp)
        E0, Ec = 0.001 * e0, 0.001 * energy(sh)
        # This implements the multi-element weighting on page 367 of 1994
        Φ0_v = Φ0(XPhi, mat, E0, Ec)
        M = sum(keys(mat)) do elm
            nonneg(mat, elm) * z(elm) / a(elm, mat) 
        end
        Φm_v = sum(keys(mat)) do elm
            ρzx_v, j, aa = ρzx(XPhi, E0, Ec, elm), 0.001 * J(Berger1982, elm), a(elm, mat)
            nonneg(mat, elm) * z(elm) / aa * Φm(XPhi, elm, j, aa, sh, E0,  ρzm(XPhi, ρzx_v, E0, Ec, elm), ρzx_v)
        end / M
        ρzx_v = sum(keys(mat)) do elm 
            nonneg(mat, elm) * z(elm) / a(elm, mat) * ρzx(XPhi, E0, Ec, elm) 
        end / M
        ρzm_v = ρzm(XPhi, ρzx_v, E0, Ec, mat)
        α_v, β_v = 0.46598 * (ρzx_v - ρzm_v), ρzm_v / sqrt(log(Φm_v / Φ0_v))
        @assert isapprox(Φm_v * exp(-(ρzm_v / β_v)^2), Φ0_v, rtol = 1.0e-6)
        @assert isapprox(exp(-((ρzx_v - ρzm_v) / α_v)^2), 0.01, atol = 1.0e-5)
        return new(sh, mat, e0, Φm_v, ρzm_v, α_v, β_v)
    end
end

"""
    ρzx(::Type{XPhi}, E0::AbstractFloat, Ec::AbstractFloat, elm::Element)

X-ray range in g/cm² from Eqn 6 Merlet 1994 or Eqn 7 in Merlet 1995
E0, Ec in keV, Z is the atomic number
"""
function ρzx(::Type{XPhi}, E0::AbstractFloat, Ec::AbstractFloat, elm::Element)::Float64
    @assert E0 > Ec "The beam energy must be larger than the edge energy!"
    @assert E0 > Ec && E0 <= 100.0
    a = 1.845e-6 * evalpoly(E0, (2.6, -0.216, 0.015, -0.000387, 0.00000501)) # Ok
    x =  2.2 - 0.0166 * E0 # Ok
    return (a * (E0^x - Ec^x) * (1.0 + 2.0 / (E0^2))) /
           ((1.078 - 0.015 * z(elm)^0.7)^(1.2 - 0.04 * E0)) # Ok
end

"""
    ρzm(::Type{XPhi}, ρzx::AbstractFloat, E0::AbstractFloat, Ec::AbstractFloat, em::Union{Element,Material})

Depth of the peak of the ϕ(ρz) curve in g/cm² (Merlet 1994 eqn 5 or Merlet 1995 eqn 7).
E0, Ec in keV, Z is the atomic number
"""
function ρzm(::Type{XPhi}, ρzx::AbstractFloat, E0::AbstractFloat, Ec::AbstractFloat, em::Union{Element,Material})::Float64
    @assert E0 < 100.0
    @assert Ec <= E0
    U0, Z = E0 / Ec, z(em)
    # return ρzx * (0.1 + 0.35 * exp(-0.07 * Z)) / (1.0 + 10.0 / (U0^10.0)) # 1994 eqn 5 
    return ρzx * (0.1 + 0.35 * exp(-0.07 * Z)) / (1.0 + 1000.0 / (Z * U0^10.0)) # 1995 eqn 7
end


"""
    τ(::Type{XPhi}, elm::Element, t::AbstractFloat)

Transmission coefficient of Zeller and Ruste 1976 (from Merlet 1994 & 1995)
"""
function τ(::Type{XPhi}, elm::Element, ρzm::AbstractFloat, ρzx::AbstractFloat)::Float64
    Z, t = z(elm), ρzm / ρzx
    τ1, τ2 = (1.0 - t)^(4.65 + 0.0356 * Z), (1.0 - t)^(1.112 + 0.00414 * Z^2)
    return τ1 + 4.65 * (τ2 - τ1) / evalpoly(Z, (3.54, 0.0356, -0.00414)) # Ok
end

mexp(::Type{XPhi}, sh::AtomicSubShell) = n(sh) == 1 ? 0.95 : 0.8

"""
    Φm(::Type{XPhi}, ::Type{XPhi}, ::Type{XPhi}, elm::Element, J::AbstractFloat, A::AbstractFloat, sh::AtomicSubShell, E0::AbstractFloat, ρzm::AbstractFloat, ρzx::AbstractFloat)

The maximum amplitude of the ϕ(ρz) curve (Located at ρzm).
Z, J, A are the atomic number, ionization potential and atomic weight (mass-fraction averaged except J which is log-averaged)
"""
function Φm(
    ::Type{XPhi}, 
    elm::Element,
    J::AbstractFloat, # Mean ionization potential (keV)
    A::AbstractFloat, # Atomic weight
    sh::AtomicSubShell,
    E0::AbstractFloat, # in keV
    ρzm::AbstractFloat,
    ρzx::AbstractFloat,
)::Float64
    @assert E0 < 100.0
    @assert J < 10.0
    Z, m, U0 = z(elm), mexp(XPhi, sh), E0 / (0.001 * energy(sh)) # Ok
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
    return τ(XPhi, elm, ρzm, ρzx) * (QUd / QU0) +
           a * (0.28 * (1.0 - 0.5 * exp(-0.1 * Z)) * b1 + 0.165 * Z^0.6 * b2) # Ok
end

"""
    Φ0(::Type{XPhi}, mat::Material, E0::AbstractFloat, Ec::AbstractFloat)

The value of the ϕ(ρz) curve at ρz=0.
"""
function Φ0(::Type{XPhi}, mat::Material, E0::AbstractFloat, Ec::AbstractFloat)::Float64
    @assert E0 < 100.0
    @assert Ec < 100.0
    Z, U0 = z(NaiveZ, mat), E0 / Ec # Merlet calls for the mass fraction averaged Z (Merlet1994, pg 367)
    @assert U0 > 1.0
    b = map( (0.02 * Z, 0.1 * Z, 0.4 * Z) ) do di
        (1.0 / di) * (1.0 - (1.0 - 1.0 / U0^di) / (di * log(U0)))
    end
    pZ = -0.25 + 0.0787 * Z^0.3 # Ok
    a = (1.87 * Z) *
        evalpoly(log(Z + 1.0), (0.0, -0.00391, 0.00721, -0.001067)) *
        (E0 / 30.0)^pZ # Ok
    return 1.0 + a * (0.27 * b[1] + (1.1 + 5.0 / Z) * (b[2] - 1.1 * b[3])) # Ok 
end


Base.range(xp::XPhi) = xp.α * sqrt(-log(0.01)) + xp.ρzm
Base.max(xp::XPhi) = xp.Φm
ϕ0(xp::XPhi) = xp.Φm * exp(-(xp.ρzm / xp.β)^2)

"""
    Φ(xphi::XPhi, ρz::AbstractFloat) # W/o absorption
    Φ(xphi::XPhi, ρz::AbstractFloat, χ::AbstractFloat) # With absorption

The antsatz for ϕ(ρz) curve in the XPhi matrix correction algorithm.
"""
function Φ(xphi::XPhi, ρz::AbstractFloat)::Float64
    ρz >= 0.0 ? #
    (
        ρz < xphi.ρzm ? #
        xphi.Φm * exp(-((ρz - xphi.ρzm) / xphi.β)^2) : #
        xphi.Φm * exp(-((ρz - xphi.ρzm) / xphi.α)^2)
    ) : #
    0.0
end
Φ(xphi::XPhi, χ::AbstractFloat, ρz::AbstractFloat)::Float64 = Φ(xphi, ρz) * exp(-ρz * χ)

"""
    ℱ(xphi::XPhi, 𝒜::AbstractFloat, ci::AbstractFloat, χ::AbstractFloat)

𝒜 - Instrumental, physical and other poorly known parameters
ci - Mass fraction of the i-th element

The generated intensity in the XPhi model.
"""
ℱ(xphi::XPhi, 𝒜::AbstractFloat, ci::AbstractFloat)::Float64 = #
    0.5 * sqrt(π) * 𝒜 * ci * xphi.Φm * (xphi.α + xphi.β * erf(xphi.ρzm / xphi.β))

"""
    ℱχp(xphi::XPhi, 𝒜::AbstractFloat, ci::AbstractFloat, χ::AbstractFloat, ρz0::AbstractFloat, ρz1::AbstractFloat)

𝒜 - Instrumental, physical and other poorly known parameters
ci - Mass fraction of the i-th element
χ - reduced mass absorption coefficient
ρz0, ρz1 - The mass depth range to integrate over

The emitted intensity in the XPhi model for the depth between ρz0 and ρz1
"""
function ℱχp(
    xphi::XPhi,
    𝒜::AbstractFloat,
    ci::AbstractFloat,
    χ::AbstractFloat,
    ρz0::AbstractFloat,
    ρz1::AbstractFloat,
)::Float64
    @assert ρz0 >= 0.0
    @assert ρz0 < ρz1
    if ρz0 < xphi.ρzm && ρz1 > xphi.ρzm
        return ℱχp(xphi, 𝒜, ci, χ, ρz0, xphi.ρzm) + ℱχp(xphi, 𝒜, ci, χ, xphi.ρzm, ρz1)
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

matrixcorrection(::Type{XPhi}, mat::Material, ashell::AtomicSubShell, e0::AbstractFloat) =
    XPhi(mat, ashell, e0)

Base.range(::Type{XPhi}, mat::Material, e0::AbstractFloat, inclDensity=false) = # 
    range(Kanaya1972, mat, e0, inclDensity)

ℱ(xphi::XPhi) = ℱ(xphi, 1.0, 1.0)

function ℱχ(xphi::XPhi, χ::AbstractFloat)::Float64
    return (0.5 * sqrt(π) * xphi.Φm) * (
        exp(χ * (0.25 * xphi.β^2 * χ - xphi.ρzm)) *
        xphi.β *
        (erf(0.5 * xphi.β * χ) + erf(xphi.ρzm / xphi.β - 0.5 * xphi.β * χ)) +
        exp(χ * (0.25 * xphi.α^2 * χ - xphi.ρzm)) * xphi.α * erfc(0.5 * xphi.α * χ)
    )
end

ℱχp(xphi::XPhi, χ::AbstractFloat, t::AbstractFloat)::Float64 =
    ℱχp(xphi, 1.0, 1.0, χ, 0.0, t)
ℱχp(xphi::XPhi, χ::AbstractFloat, t0::AbstractFloat, t1::AbstractFloat)::Float64 =
    ℱχp(xphi, 1.0, 1.0, χ(material(xphi), xray, θtoa), t0, t1)

ϕ(xphi::XPhi, ρz) = Φ(xphi, ρz)
