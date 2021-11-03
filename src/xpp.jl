using PeriodicTable
using DataFrames
# Implements Pouchou and Pichoir's XPP model from the description in the book
# Electron Probe Quantification edited by Kurt Heinrich and Dale E. Newbury

struct XPP <: MatrixCorrection
    # Configuration data
    subshell::AtomicSubShell
    material::Material
    E0 # Beam energy (eV)
    # Computed items
    ϕ0 # ϕ(ρz) at the surface
    A # Amplitude  factor
    a # Width factor
    B # Amplitude  factor
    b # Width factor
    F # Integral under the ϕ(ρz) curve

    """
        XPP(mat::Material, ashell::AtomicSubShell, E0::AbstractFloat)

    Construct an XPP object for the specified material, atomicsubshell,
    beam energy (in eV).
    """
    function XPP(mat::Material, ashell::AtomicSubShell, E0::AbstractFloat)
        # XPP calculations expect E0, eLv, J in keV
        e0, Mv, Jv, = 0.001 * E0, M(XPP, mat), 0.001 * J(XPP, mat)
        Zbarbv, eLv = Zbarb(XPP, mat), 0.001 * energy(ashell)
        U0v, V0v, ηbarv = e0 / eLv, e0 / Jv, ηbar(XPP, Zbarbv)
        @assert U0v > 1.0 "The beam energy must be larger than the subshell edge energy."
        @assert V0v > 1.0 "The beam energy must be larger than the mean energy loss. ($(mat), $(ashell), $(E0))"
        Dv, Pjv, mv, Wbarv = D(XPP, Jv), P(XPP, Jv), m(XPP, ashell), Wbar(XPP, ηbarv)
        invSv = invS(XPP, U0v, V0v, Mv, Dv, Pjv, T(XPP, Pjv, mv))
        qv, Rv, ϕ0v, QlAv = q(XPP, Wbarv), R(XPP, ηbarv, Wbarv, U0v), ϕ0(XPP, U0v, ηbarv), QlA(XPP, U0v, eLv, mv)
        Fv, FoverRbarv = F(XPP, Rv, invSv, QlAv), FoverRbar(XPP, X(XPP, Zbarbv), Y(XPP, Zbarbv), U0v)
        Rbarv = Rbar(XPP, Fv, FoverRbarv, ϕ0v)
        bv = b(XPP, Rbarv, ϕ0v, Fv)
        Pv = P(XPP, Zbarbv, U0v, Fv, Rbarv, ϕ0v, bv)
        ϵv = ϵ(XPP, Pv, bv, ϕ0v, Fv, Rbarv)
        av, Bv = a(XPP, bv, ϵv), B(XPP, bv, Fv, Pv, ϕ0v, ϵv)
        Av = A(XPP, Bv, bv, ϕ0v, Fv, ϵv)
        return new(ashell, mat, E0, ϕ0v, Av, av, Bv, bv, Fv)
    end
end

"""
    M(mat::Material)

P&P's M from top line in page 35
"""
M(::Type{XPP}, mat::Material) = #C1
    sum(NeXLCore.nonneg(mat, elm) * z(elm) / a(elm, mat) for elm in keys(mat))

"""
    J(::Type{XPP}, mat::Material)

Mean ionization potential for the specified material in eV. (PAP1991 Eqn 6)
"""
NeXLCore.J(::Type{XPP}, mat::Material) = #C1
    exp(sum(NeXLCore.nonneg(mat, elm) * (z(elm) / a(elm, mat)) * log(J(Zeller1973, elm)) for elm in keys(mat)) / M(XPP, mat))


# From PAP eqn 8
D(::Type{XPP}, J) = (6.6e-6, 1.12e-5 * (1.35 - 0.45 * J^2), 2.2e-6 / J) #C1

# From PAP eqn 8
P(::Type{XPP}, J) = (0.78, 0.1, -0.5 + 0.25 * J) #C1

# From PAP near Eqn 11
T(::Type{XPP}, P, m) = collect(1.0 + P[k] - m for k = 1:3) #C1

"""
    dEdρs(args::Dict{Label,AbstractFloat}, mat::MaterialLabel, Ekev::AbstractFloat, elms)

The function P&P uses to describe the deceleration of an electron in a material.
Output units are (keV/cm)/(g/cm^3) = keV cm^2/g. (PAP eqn 5)
"""
function dEdρs(::Type{XPP}, mat::Material, Ekev::AbstractFloat) #C1
    @assert Ekev < 50.0 "It appears that the beam energy is in keV not eV. ($Ekev)"
    Jkev = 0.001 * J(XPP, mat) # XPP expects in keV
    # PAP Eqn 8
    function f(mat, Ekev, Jkev)
        d, p, v = D(XPP, j), P(XPP, j), Ekev / Jkev
        return sum(i -> d[i] * v^p[i], 1:3)
    end
    return -(M(XPP, mat) / Jkev) / f(mat, Ekev, Jkev)
end

"""
    R0(J, D, P, M, Ekev)

Total trajectory (range) of an electron with initial energy Ekev. (in cm/(g/cm^3))
"""
R0(::Type{XPP}, J, D, P, M, Ekev) = #C1
    sum(J^(1.0 - P[k]) * D[k] * Ekev^(1.0 + P[k]) / (1.0 + P[k]) for k = 1:3) / M

"""
    ϕxpp(ρz, A, a, B, b, ϕ0)

Compute the shape of the ϕ(ρz) curve in the XPP model.
"""
ϕxpp(::Type{XPP}, ρz, A, a, B, b, ϕ0) = #C1
    A * exp(-a * ρz) + (B * ρz + ϕ0 - A) * exp(-b * ρz)

"""
    invS(U0, V0, M, D, P, T)

Computes 1/S where S is the stopping power.
"""
invS(::Type{XPP}, U0, V0, M, D, P, T) = #C1
    U0 / (V0 * M) * sum(D[k] * ((V0 / U0)^P[k]) * (T[k] * U0^T[k] * log(U0) - U0^T[k] + 1.0) / (T[k]^2) for k = 1:3)

"""
    S(mat, ashell, E)

Computes S, the stopping power at the specified energy (in keV)
"""
function S(::Type{XPP}, mat::Material, ashell::AtomicSubShell, Ekev)
    @assert Ekev < 50.0 "It appears that the beam energy is in keV not eV. ($Ekev)"
    jkev = 0.001 * J(XPP, mat) # XPP expects in keV
    return 1.0 / invS(XPP, Ekev / (0.001 * energy(ashell)), Ekev / jkev, M(XPP, mat), D(XPP, jkev), P(XPP, jkev), T(XPP, P(XPP, jkev), m(XPP, ashell)))
end

"""
    QlA(U,El,m)

Computes the relative ionization cross-section.
"""
QlA(::Type{XPP}, U, El, m) = log(U) / ((U^m) * (El^2)) #C1


"""
    m(ashell::AtomicSubShell)

Returns the ionization cross-section exponent for QlA(U,El,m(ashell))
"""
function m(::Type{XPP}, ashell::AtomicSubShell) #C1
    if isequal(n(ashell), 1)
        return 0.86 + 0.12 * exp(-(ashell.z / 5.0)^2)
    else
        return n(ashell) == 2 ? 0.82 : 0.78
    end
end

# PAP appendix 1
Zbarb(::Type{XPP}, mat::Material) = #C1
    sum(NeXLCore.nonneg(mat, elm) * sqrt(z(elm)) for elm in keys(mat))^2

# PAP appendix 1
ηbar(::Type{XPP}, Zbarb) = #C2
    1.75e-3 * Zbarb + 0.37 * (1.0 - exp(-0.015 * Zbarb^1.3))

# PAP appendix 1
Wbar(::Type{XPP}, ηbar) = #C2
    0.595 + ηbar / 3.7 + ηbar^4.55

q(::Type{XPP}, Wbar) = #C1
    (2.0 * Wbar - 1.0) / (1.0 - Wbar)

# PAP appendix 1
JU0(::Type{XPP}, U0) = #C2
    1.0 + U0 * (log(U0) - 1.0)

# PAP appendix 1
G(::Type{XPP}, U0, q) = #C2
    (U0 - 1.0 - (1.0 - 1.0 / U0^(1.0 + q)) / (1.0 + q)) / ((2.0 + q) * JU0(XPP, U0))

# PAP appendix 1
R(::Type{XPP}, ηbar, Wbar, U0) = #C2
    1.0 - ηbar * Wbar * (1.0 - G(XPP, U0, q(XPP, Wbar)))

"""
    R(mat::Material, u0)

Backscatter factor as a function of material and overvoltage.
Compares well to PAP Figure 23.
"""
function R(::Type{XPP}, mat::Material, u0)
    zbarbv = Zbarb(XPP, mat)
    ηbarv = ηbar(XPP, zbarbv)
    wbar = Wbar(XPP, ηbarv)
    return R(XPP, ηbarv, wbar, u0)
end

"""
   ϕ0(U0, ηbar)

The value of the ϕ(ρz) curve at ρz=0.0.
"""
ϕ0(::Type{XPP}, U0, ηbar) = #C1
    1.0 + 3.3 * (1.0 - 1.0 / (U0^(2.0 - 2.3 * ηbar))) * ηbar^1.2

X(::Type{XPP}, Zbarb) = #C3
    1.0 + 1.3 * log(Zbarb)

Y(::Type{XPP}, Zbarb) = #C3
    0.2 + Zbarb / 200.0

FoverRbar(::Type{XPP}, X, Y, U0) = #C3
    1.0 + X * log(1.0 + Y * (1.0 - 1.0 / (U0^0.42))) / log(1.0 + Y)

# PAP eqn 13
F(::Type{XPP}, R, invS, QlA) = #C1
    R * invS / QlA

# PAP eqn 28
Rbar(::Type{XPP}, F, FoverRbar, ϕ0) = # C1
    FoverRbar >= ϕ0 ? F / FoverRbar : F / ϕ0

# PAP eqn 28
g(::Type{XPP}, Zbarb, U0) = #C1
    0.22 * log(4.0 * Zbarb) * (1.0 - 2.0 * exp(-Zbarb * (U0 - 1.0) / 15.0))

# PAP eqn 28
h(::Type{XPP}, Zbarb, U0) = #C1
    1.0 - 10.0 * (1.0 - 1.0 / (1.0 + U0 / 10.0)) / (Zbarb^2)

"""
    b(Rbar, ϕ0, F)

XPP ϕ(ρz) model parameter b
""" # PAP Appendix 4
b(::Type{XPP}, Rbar, ϕ0, F) = #C1
    sqrt(2.0) * (1.0 + sqrt(1.0 - Rbar * ϕ0 / F)) / Rbar

# PAP eqn 29
P(::Type{XPP}, Zbarb, U0, F, Rbar, ϕ0, b) = #C1
    min(g(XPP, Zbarb, U0) * h(XPP, Zbarb, U0)^4, 0.9 * b * Rbar^2 * (b - 2.0 * ϕ0 / F)) * F / (Rbar^2)

function ϵ(::Type{XPP}, P, b, ϕ0, F, Rbar) #C1
    a = (P + b * (2.0 * ϕ0 - b * F)) / (b * F * (2.0 - b * Rbar) - ϕ0)
    tmp = (a - b) / b
    # print("a = ",a,"(a-b)/b = ",tmp)
    return abs(tmp) > 1.0e-6 ? tmp : sign(tmp) * 1.0e-6
end

"""
    a(b, ϵ)

XPP ϕ(ρz) model parameter a
"""
NeXLCore.a(::Type{XPP}, b, ϵ) = #C1
    b * (1.0 + ϵ)

"""
    B(b, F, P, ϕ0, ϵ)

XPP ϕ(ρz) model parameter B
"""
B(::Type{XPP}, b, F, P, ϕ0, ϵ) = #C1
    (b^2 * F * (1.0 + ϵ) - P - ϕ0 * b * (2.0 + ϵ)) / ϵ

"""
    A(B, b, ϕ0, F, ϵ)

XPP ϕ(ρz) model parameter A
"""
A(::Type{XPP}, B, b, ϕ0, F, ϵ) = #C1
    (B / b + ϕ0 - b * F) * (1 + ϵ) / ϵ

function Fχ(::Type{XPP}, χ, A, a, B, b, ϕ0)
    ϵ = (a - b) / b
    return (ϕ0 + B / (b + χ) - A * b * ϵ / (b * (1.0 + ϵ) + χ)) / (b + χ)
end


"""
    Fχp(::Type{XPP}, χ, A, a, B, b, ϕ0, τ)

The integral of the ϕ(ρz) exp(-χ ρz) from 0 to τ.
"""
Fχp(::Type{XPP}, χ, A, a, B, b, ϕ0, τ) =
    (A * (1.0 - exp(-(τ * (a + χ))))) / (a + χ) +
    (A * (-1.0 + exp(-(τ * (b + χ))))) / (b + χ) +
    ((1 - exp(-(τ * (b + χ)))) * ϕ0) / (b + χ) +
    (B * (-1 + exp(τ * (b + χ)) - τ * (b + χ))) / (exp(τ * (b + χ)) * ((b + χ)^2))

"""
Represents the essential intermediary values for an XPP matrix correction of
characteristic X-rays from a particular atomic sub-shell in a particular material.
"""
function NeXLUncertainties.asa(::Type{DataFrame}, xpps::XPP...)
    return DataFrame(
        SubShell = [xpp.subshell for xpp in xpps],
        Material = [name(xpp.material) for xpp in xpps],
        BeamEnergy = [xpp.E0 for xpp in xpps],
        Phi0 = [xpp.ϕ0 for xpp in xpps],
        A = [xpp.A for xpp in xpps],
        a = [xpp.a for xpp in xpps],
        B = [xpp.B for xpp in xpps],
        b = [xpp.b for xpp in xpps],
        F = [xpp.F for xpp in xpps],
    )
end

Base.show(io::IO, xpp::XPP) = print(io, "XPP[$(xpp.subshell) in $(name(xpp.material)) at $(0.001*xpp.E0) keV]")

Fχ(xpp::XPP, xray::CharXRay, θtoa::Real) = Fχ(XPP, χ(material(xpp), xray, θtoa), xpp.A, xpp.a, xpp.B, xpp.b, xpp.ϕ0)


Fχp(xpp::XPP, xray::CharXRay, θtoa::Real, τ::Real) =
    Fχp(XPP, χ(material(xpp), xray, θtoa), xpp.A, xpp.a, xpp.B, xpp.b, xpp.ϕ0, τ)

F(xpp::XPP) = xpp.F

NeXLCore.minproperties(::Type{XPP}) = (:BeamEnergy, :TakeOffAngle)

"""
    ϕ(ρz, xpp::XPP)

Computes the ϕ(ρz) curve according to the XPP algorithm.
"""
ϕ(xpp::XPP, ρz) = ϕxpp(XPP, ρz, xpp.A, xpp.a, xpp.B, xpp.b, xpp.ϕ0)


"""
    range(::Type{XPP}, mat::MaterialLabel, e0, inclDensity=true)
    range(ty::Type{<:BetheEnergyLoss}, mat::Material, e0::Float64, inclDensity=true; emin=50.0, mip::Type{<:NeXLMeanIonizationPotential}=Berger1982)

Total trajectory (range) of an electron with initial energy e0 (eV). (Units: inclDensity ? cm : cm/(g/cm³))
"""
function Base.range(::Type{XPP}, mat::Material, e0::Real, inclDensity=true)
    j = 0.001 * J(XPP, mat) # XPP expects in keV
    return R0(XPP, j, D(XPP, j), P(XPP, j), M(XPP, mat), 0.001 * e0) / (inclDensity ? density(mat) : 1.0)
end

"""
    matrixcorrection(::Type{XPP}, mat::Material, ashell::AtomicSubShell,e0)

Constructs an instance of the XPP algorithm.
"""
matrixcorrection(::Type{XPP}, mat::Material, ashell::AtomicSubShell, e0) = XPP(mat, ashell, e0)
