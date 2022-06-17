# Generic methods for matrix correctio

using DataFrames
using PeriodicTable
using Roots

"""
`MatrixCorrection` is an abstract type for computing ϕ(ρz)-type matrix correction algorithms.  A sub-class
`MCA <: MatrixCorrection` should implement

    # Integral of the area under the ϕ(ρz)-curve
    ℱ(mc::MCA)
    # Integral of the transmitted area under the ϕ(ρz)-curve
    ℱχ(mc::MCA, xray::CharXRay, θtoa::AbstractFloat)
    # Area under 0 to τ under the transmitted ϕ(ρz)-curve
    ℱχp(mc::MCA, xray::CharXRay, θtoa::AbstractFloat, τ::AbstractFloat)
    # The ϕ(ρz)-curve
    ϕ(mc::MCA, ρz)
    # A factory method from MCA
    matrixcorrection(::Type{MCA}, mat::Material, ashell::AtomicSubShell, e0)::MCA

The algorithm should precompute as much as possible based on the `Material`, `AtomicSubShell` and beam
energy.  The class `MCA` should define member variables `subshell::AtomicSubShell`, `material::Material` and
`E0::AbstractFloat` to store the input values. See `XPP` for an example.

From these methods, other methods like `Z(...)`, `A(...)` are implemented.
"""
abstract type MatrixCorrection end

"""
    ℱ(mc::MatrixCorrection)

Integral of the ϕ(ρz)-curve from ρz = 0 to ∞.
"""
ℱ(mc::MatrixCorrection) = error("$mc does not implement ℱ(mc::$mc)")
"""
    ℱχ(mc::MatrixCorrection, χ::AbstractFloat)
    ℱχ(mc::MatrixCorrection, xray::CharXRay, θtoa::AbstractFloat)

Integral of the area under the absorption corrected ϕ(ρz)-curve from ρz = 0 to ∞.
"""
ℱχ(mc::MatrixCorrection, ::AbstractFloat)  = error("$mc does not implement ℱχ(mc::$mc, χ::AbstractFloat)")
ℱχ(mc::MatrixCorrection, xray::CharXRay, θtoa::AbstractFloat)  = ℱχ(mc, χ(material(mc), xray, θtoa))

"""
    ℱχp(mc::NeXLMatrixCorrection, χ::AbstractFloat, τ::AbstractFloat)
    ℱχp(mc::NeXLMatrixCorrection, xray::CharXRay, θtoa::AbstractFloat, τ::AbstractFloat)
    ℱχp(mc::MatrixCorrection, xray::CharXRay, θtoa::AbstractFloat, t0::AbstractFloat, t1::AbstractFloat)

The partial integral of the absorption corrected ϕ(ρz) curve from ρz = 0 to τ or from t0 to t1
"""
ℱχp(mc::MatrixCorrection, ::AbstractFloat, ::AbstractFloat) = #
    error("$mc does not implement ℱχp(mc::$mc, χ::AbstractFloat, τ::AbstractFloat)")
ℱχp(mc::MatrixCorrection, cxr::CharXRay, θtoa::AbstractFloat, t::AbstractFloat)  = #
    ℱχp(mc, χ(material(mc), cxr, θtoa), t)
ℱχp(mc::MatrixCorrection, xray::CharXRay, θtoa::AbstractFloat, t0::AbstractFloat, t1::AbstractFloat) = #
    ℱχp(mc, xray, θtoa, t1) - ℱχp(mc, xray, θtoa, t0)

"""
    ϕ(mc::MatrixCorrection, ρz)

The ϕ(ρz)-curve for visualization and other purposes.
"""
ϕ(mc::MatrixCorrection, ρz)  = error("$mc does not implement ϕ(mc::$mc, ρz)")

"""
    NeXLCore.atomicsubshell(mc::MatrixCorrection)

The sub-shell for which this `MatrixCorrection` has been calculated.
"""
NeXLCore.atomicsubshell(mc::MatrixCorrection) = mc.subshell

"""
    NeXLCore.material(mc::MatrixCorrection) = mc.material

The material for which this `MatrixCorrection` has been calculated.
"""
NeXLCore.material(mc::MatrixCorrection) = mc.material

"""
    beamEnergy(mc::MatrixCorrection) = mc.E0

The beam energy (eV) for which this `MatrixCorrection` has been calculated.
"""
beamEnergy(mc::MatrixCorrection) = mc.E0


NeXLCore.minproperties(::Type{<:MatrixCorrection}) = (:BeamEnergy, :TakeOffAngle)

"""
    χ(mat::Material, xray::CharXRay, θtoa)

Angle adjusted mass absorption coefficient.
"""
χ(mat::Material, xray::CharXRay, θtoa::AbstractFloat) = mac(mat, xray) * csc(θtoa)
χ(mat::Material, ea::AbstractFloat, θtoa::AbstractFloat) = mac(mat, ea) * csc(θtoa)

"""
    ϕabs(mc::MatrixCorection, ρz, xray::CharXRay, θtoa::AbstractFloat)
    ϕabs(mc::MatrixCorrection, ρz::AbstractFloat, χ::AbstractFloat)

Computes the absorbed ϕ(ρz) curve according to the XPP algorithm.
"""
ϕabs(mc::MatrixCorrection, ρz::AbstractFloat, xray::CharXRay, θtoa::AbstractFloat) = ϕabs(mc, χ(material(mc), xray, θtoa), ρz)
ϕabs(mc::MatrixCorrection, ρz::AbstractFloat, χ::AbstractFloat) = ϕ(mc, ρz) * exp(-χ * ρz)

"""
    ZA(
      unk::MatrixCorrection,
      std::MatrixCorrection,
      xray::CharXRay,
      θunk::AbstractFloat,
      θstd::AbstractFloat
    )

The atomic number (Z) and absorption (A) correction factors.
"""
function ZA(unk::MatrixCorrection, std::MatrixCorrection, xray::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat)
    @assert isequal(unk.subshell, inner(xray)) "Unknown and X-ray don't match"
    @assert isequal(std.subshell, inner(xray)) "Standard and X-ray don't match"
    return ℱχ(unk, xray, θunk) / ℱχ(std, xray, θstd)
end

"""
    Z(unk::MatrixCorrection, std::MatrixCorrection)

The atomic number correction factor.
"""
function Z(unk::MatrixCorrection, std::MatrixCorrection)
    @assert isequal(unk.subshell, std.subshell)
    "Unknown and standard matrix corrections don't apply to the same sub-shell."
    return ℱ(unk) / ℱ(std)
end


"""
    A(unk::MatrixCorrection, std::MatrixCorrection, xray::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat)
    A(unk::MatrixCorrection, xray::CharXRay, θ::AbstractFloat)

The absorption correction factors.
"""
function A(unk::MatrixCorrection, std::MatrixCorrection, xray::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat)
    @assert isequal(unk.subshell, inner(xray)) "Unknown and X-ray don't match"
    @assert isequal(std.subshell, inner(xray)) "Standard and X-ray don't match"
    return ZA(unk, std, xray, θunk, θstd) / Z(unk, std)
end

function A(unk::MatrixCorrection, xray::CharXRay, θunk::AbstractFloat)
    return ℱχ(unk, xray, θunk)/ℱ(unk)
end

function correctcontinuum(mc::MatrixCorrection, θtoa::AbstractFloat)
    @assert isnothing(atomicsubshell(mc))  "Use only for continuum correction"
    return ℱχ(mc, χ(material(mc), edgeenergy(mc), θtoa)) / ℱ(mc)
end

"""
    kcoating(ty::Type{<:MatrixCorrection}, subtrate::Material, coating::Material, cxr::CharXRay, e0::AbstractFloat, toa::AbstractFloat, τ::AbstractFloat)

Estimate the k-ratio for a coating of mass-thickness τ (g/cm²) on the specified substrate. The standard for the coating
is assumed to be of the same material as the coating.
"""
function kcoating(
    ty::Type{<:MatrixCorrection},
    subtrate::Material,
    coating::Material,
    cxr::CharXRay,
    e0::AbstractFloat,
    toa::AbstractFloat,
    τ::AbstractFloat,
)
    @assert element(cxr) in keys(coating) "$cxr must be produced by one of the elements in $(name(coating))"
    @assert e0 > energy(inner(cxr)) "The beam energy must exceed the edge energy for $cxr."
    submc, coatmc = matrixcorrection(ty, subtrate, inner(cxr), e0), matrixcorrection(ty, coating, inner(cxr), e0)
    return ℱχp(submc, cxr, toa, τ) / ℱχ(coatmc, cxr, toa)
end

""""
    massthickness(ty::Type{<:MatrixCorrection}, subtrate::Material, coating::Material, cxr::CharXRay, e0::AbstractFloat, toa::AbstractFloat, k::AbstractFloat)

Estimate the mass-thickness of a ultra-thin layer of a `coating` material on a `substrate` from a measured k-ratio `k`
of a characteristic X-ray `cxr`.  Works for k-ratios of the order of 1 %.  The standard for the coating is assumed
to be of the same material as the coating.
"""
function NeXLCore.massthickness(
    ty::Type{<:MatrixCorrection},
    substrate::Material,
    coating::Material,
    cxr::CharXRay,
    e0::AbstractFloat,
    toa::AbstractFloat,
    k::AbstractFloat,
)
    submc, coatingmc = matrixcorrection(ty, substrate, inner(cxr), e0), matrixcorrection(ty, coating, inner(cxr), e0)
    f(τ) = k - ℱχp(submc, cxr, toa, τ) / ℱχ(coatingmc, cxr, toa)
    return find_zero(f, 0.01 * range(ty, substrate, e0, false), Roots.Order1())
end

function coatingasfilm(
    ty::Type{<:MatrixCorrection},
    substrate::Material,
    coating::Material,
    cxr::CharXRay,
    e0::AbstractFloat,
    toa::AbstractFloat,
    k::AbstractFloat,)
    @assert haskey(coating.properties, :Density)
    @assert element(cxr) in keys(coating) "The element associated with $cxr is not present in $coating."
    return Film(coating, massthickness(ty, substrate, coating, cxr, e0, toa, k) / coating[:Density])
end


"""
    correctkratios(krs::AbstractVector{KRatio}, coating::Material, θtoa::AbstractFloat, ρz::AbstractFloat)::Vector{KRatio}

This function is mainly for pedagogical purposes.  It takes a `KRatio[]`, a coating `Material` on the unknown,
and a mass-thickness (g/cm²) and creates a new `KRatio[]` that accounts for the intensity missing due to absorption
by the coating.  Favor the coating correction built into `ZAFCorrection` or `MultiZAF`.
"""
function correctkratios(krs::AbstractVector{KRatio}, coating::Material, ρz::AbstractFloat)::Vector{KRatio}
    res = KRatio[]
    for kr in krs
        if !(kr.element in keys(coating))
            kratio = kr.kratio / exp(-χ(coating, brightest(kr.xrays), kr.unkProps[:TakeOffAngle]) * ρz)
            push!(res, KRatio(kr.xrays, kr.unkProps, kr.stdProps, kr.standard, kratio))
        else
            push!(res, kr)
        end
    end
    return res
end
