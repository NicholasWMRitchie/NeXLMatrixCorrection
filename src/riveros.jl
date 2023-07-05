using SpecialFunctions: erfc, erf

"""
@article{riveros1993review,
  title={Review of ϕ(ρz) curves in electron probe microanalysis},
  author={Riveros, Jose and Castellano, Gustavo},
  journal={X-Ray Spectrometry},
  volume={22},
  number={1},
  pages={3--10},
  year={1993},
  publisher={Wiley Online Library}
}

The instruction in Packwood1991 in Electron Probe Quantitation is:
\"For compounds weight averaging is used for all appropriate variables: Z, Z/A, η, and (Z/A)log(1.166(E0-Ec)/2J)\"
"""
struct Riveros1993 <: MatrixCorrection
    material::Material{Float64, Float64}
    E0::Float64  # Beam energy
    subshell::Union{Nothing,AtomicSubShell}
    Ea::Float64 # Sub-shell energy
    γ::Float64
    ϕ0::Float64
    α::Float64
    β::Float64

    function Riveros1993(matp::Material, ashell::Union{AtomicSubShell,Nothing}, ea::AbstractFloat, e0::AbstractFloat)
        mat = convert(Material{Float64,Float64}, matp)
        e0k, eck, u0 = 0.001 * e0, 0.001 * ea, e0 / ea
        @assert u0 >= 1.0
        αz(elm) = (2.14e5*z(elm)^1.16/(a(elm)*e0k^1.25))*sqrt(log(1.166*e0/J(Berger1982,elm))/(e0k-eck))
        βz(elm) = (1.1e5*z(elm)^1.5)/((e0k-eck)*a(elm))
        ηm = zfractionaverage(el->η(el, e0), mat) # Use Donovan's averaging (not weight averaging...)
        ϕ0 = 1.0 + (ηm*u0*log(u0))/(u0-1.0) # ok
        γ = (1.0 + ηm)*(u0*log(u0))/(u0-1.0) # ok
        α = zfractionaverage(αz, mat)
        β = zfractionaverage(βz, mat)
        @assert α > 0.0
        @assert β > 0.0
        return new(mat, e0, ashell, ea, γ, ϕ0, α, β)
    end
end

ϕ(rv::Riveros1993, ρz) = exp(-(rv.α*ρz)^2)*(rv.γ - (rv.γ-rv.ϕ0)*exp(-rv.β*ρz)) # ok

function ℱ(rv::Riveros1993)::Float64
    gg = rv.β/(2.0*rv.α)
    return gg < 22.3 ? (√π*(rv.γ - exp(gg^2)*(rv.γ - rv.ϕ0)*erfc(gg)))/(2.0*rv.α) : 1.0
end

function ℱχ(rv::Riveros1993, χm::AbstractFloat)::Float64
    ff, gg = χm/(2.0*rv.α), (rv.β + χm)/(2.0*rv.α)
    return gg < 22.3 ? (sqrt(π)*(
            exp(ff^2)*rv.γ*erfc(ff) - exp(gg^2)*(rv.γ - rv.ϕ0)*erfc(gg)))/ #
            (2.0*rv.α) : 0.0
end

function ℱχp(rv::Riveros1993, xray::CharXRay, θtoa::AbstractFloat, τ::AbstractFloat)::Float64
    @assert isnothing(rv.subshell) || (inner(xray) == rv.subshell)
    χm = χ(material(rv), xray, θtoa)
    @assert rv.α > 0.0 && χm > 0.0
    ff, gg = χm/(2.0*rv.α), (rv.β + χm)/(2.0*rv.α)
    return gg < 22.3 ? (sqrt(π)*rv.α*(exp(ff^2)*rv.γ*
            (-((χm*erf(ff^2))/χm) +
            ((2*rv.α^2*τ + χm)*erf((2.0*rv.α^2*τ + χm)/(2.0*rv.α)))/ (2*rv.α^2*τ + χm)) -
           exp(gg^2)*rv.γ*
            (-(((rv.β + χm)*erf(gg))/abs(rv.β + χm)) +
              ((rv.β + 2.0*rv.α^2*τ + χm)*erf(abs(rv.β + 2*rv.α^2*τ + χm))/(2.0*rv.α))/
               abs(rv.β + 2.0*rv.α^2*τ + χm))))/(2.0*rv.α^2) : 0.0
end

NeXLCore.edgeenergy(rv::Riveros1993) = rv.Ea

matrixcorrection(::Type{Riveros1993}, mat::Material, ashell::AtomicSubShell, e0::AbstractFloat) = Riveros1993(mat, ashell, energy(ashell), e0)

continuumcorrection(::Type{Riveros1993}, mat::Material, ea::AbstractFloat, e0::AbstractFloat) = Riveros1993(mat, nothing, ea, e0)

Base.range(::Type{Riveros1993}, mat::Material, e0::AbstractFloat, inclDensity=true) = range(XPP, mat, e0, inclDensity)

NeXLCore.minproperties(::Type{Riveros1993}) = (:BeamEnergy, :TakeOffAngle)
