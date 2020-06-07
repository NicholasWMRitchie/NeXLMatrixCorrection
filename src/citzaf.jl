# Implements Armstrong's CitZAF algorithm as presented in the Green Book (Heinrich, 1991)

struct CitZAF <: MatrixCorrection
    material::Material
    E0::Float64  # Beam energy
    subshell::Union{Nothing,AtomicSubShell}
    Ea::Float64 # Sub-shell energy
    γ0::Float64
    q::Float64
    α::Float64
    β::Float64

    function CitZAF(mat::Material, ashell::Union{AtomicSubShell,Nothing}, ea::Float64, e0::Float64)
        function ϕ0L(u0, ηbar)
            a = 3.43378 + (-10.7872 + (10.97628 - 3.62286 / u0) / u0) / u0
            b = -0.59299 + (21.55329 + (-30.55248 + 9.59218 / u0) / u0) / u0
            return 1.0 + (ηbar / (1.0 + ηbar)) * (a + b * log(1.0 + ηbar))
        end
        function ηbarLS(z, e0k)
            h1 = 1.0e-4 * (-52.3791 + z * (150.48371 + z * (-1.67373 + z * 0.00716)))
            h2 = 1.0e-4 * (-1112.8 + z * (30.289 - z * 0.15498))
            return h1 * (1.0 + h2 * log(e0k / 20.0))
        end
        e0k, eck, u0 = 0.001 * e0, 0.001 * ea, e0 / ea
        @assert u0 >= 1.0
        zbar =
            sum(value(mat[elm]) * z(elm) / a(elm, mat) for elm in keys(mat)) /
            sum(value(mat[elm]) / a(elm, mat) for elm in keys(mat))
        abar = sum(value(mat[elm]) for elm in keys(mat)) / sum(value(mat[elm]) / a(elm, mat) for elm in keys(mat))
        γ0 = (5π * u0) / ((u0 - 1.0) * log(u0)) * (log(u0) - 5.0 + 5.0 * u0^-0.2)
        ηbar = mapreduce(elm -> mat[elm] * ηbarLS(z(elm), e0k), +, keys(mat))
        ϕ0 = ϕ0L(u0, ηbar) # 1.0 + 2.8*(1.0 -  0.9/u0)*ηbar
        q = (γ0 - ϕ0) / γ0
        α = 2.97e5 * (zbar^1.05 / (abar * e0k^1.25)) * sqrt(log(1.166 * e0 / J(XPP, mat)) / (e0k - eck))
        β = (8.5e5 * zbar^2) / (abar * e0k^2 * (γ0 - 1.0))
        return new(mat, e0, ashell, ea, γ0, q, α, β)
    end
end

ϕ(cz::CitZAF, ρz) = cz.γ0 * exp(-(cz.α * ρz)^2) * (1.0 - cz.q * exp(-cz.β * ρz))

F(cz::CitZAF) = ((cz.α - cz.q * cz.α + cz.β) * cz.γ0) / (cz.α * (cz.α + cz.β))

function Fχ(cz::CitZAF, ea::Float64, θtoa::Real)
    @assert isnothing(cz.subshell)  "Use only for continuum correction"
    χm = χ(material(cz), ea, θtoa)
    return cz.γ0 * (1.0 / (cz.α + χm) - cz.q / (cz.α + cz.β + χm))
end

function Fχ(cz::CitZAF, xray::CharXRay, θtoa::Real)
    @assert !isnothing(cz.subshell) "Use only for characteristic correction"
    @assert inner(xray) == cz.subshell
    χm = χ(material(cz), xray, θtoa)
    return cz.γ0 * (1.0 / (cz.α + χm) - cz.q / (cz.α + cz.β + χm))
end

function Fχp(cz::CitZAF, xray::CharXRay, θtoa::Real, τ::Real)
    @assert isnothing(cz.subshell) || (inner(xray) == cz.subshell)
    χm = χ(material(cz), xray, θtoa)
    return ((1.0 - exp(-τ * (cz.α + χm))) * cz.γ0) / (cz.α + χm) +
           ((-1.0 + exp(-τ * (cz.α + cz.β + χm))) * cz.q * cz.γ0) / (cz.α + cz.β + χm)
end

matrixcorrection(::Type{CitZAF}, mat::Material, ashell::AtomicSubShell, e0::Float64) = CitZAF(mat, ashell, energy(ashell), e0)

continuumcorrection(::Type{CitZAF}, mat::Material, ea::Float64, e0::Float64) = CitZAF(mat, nothing, ea, e0)

Base.range(::Type{CitZAF}, mat::Material, e0::Float64) = range(XPP, mat, e0)

NeXLCore.minproperties(::Type{CitZAF}) = (:BeamEnergy, :TakeOffAngle)
