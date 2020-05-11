# Generic methods for matrix correctio

using DataFrames
using PeriodicTable
using Roots

"""
`MatrixCorrection` is an abstract type for computing ϕ(ρz)-type matrix correction algorithms.  Subclasses should
implement

    F(mc::MatrixCorrection) # Integral of the area under the ϕ(ρz)-curve
    Fχ(mc::MatrixCorrection, xray::CharXRay, θtoa::Real) # Integral of the transmitted area under the ϕ(ρz)-curve
    Fχp(mc::NeXLMatrixCorrection, xray::CharXRay, θtoa::Real, τ::Real) # Area under 0 to τ under the transmitted ϕ(ρz)-curve
    ϕ(mc::MatrixCorrection, ρz) # The ϕ(ρz)-curve

From these methods, other methods like `Z(...)`, `A(...)` are implemented.
"""
abstract type MatrixCorrection end


NeXLCore.atomicsubshell(mc::MatrixCorrection) = mc.subshell
NeXLCore.material(mc::MatrixCorrection) = mc.material
beamEnergy(mc::MatrixCorrection) = mc.E0

"""
    χ(mat::Material, xray::CharXRay, θtoa)
Angle adjusted mass absorption coefficient.
"""
χ(mat::Material, xray::CharXRay, θtoa::AbstractFloat) = mac(mat, xray) * csc(θtoa)

"""
    ϕabs(mc::MatrixCorection, ρz, xray::CharXRay, θtoa::AbstractFloat)

Computes the absorbed ϕ(ρz) curve according to the XPP algorithm.
"""
ϕabs(mc::MatrixCorrection, ρz, xray::CharXRay, θtoa::AbstractFloat) = ϕ(mc, ρz) * exp(-χ(material(mc), xray, θtoa) * ρz)

"""
    ZA(
      unk::MatrixCorrection,
      std::MatrixCorrection,
      xray::CharXRay,
      θunk::AbstractFloat,
      θstd::AbstractFloat
    )

The atomic number and absorption correction factors.
"""
function ZA(unk::MatrixCorrection, std::MatrixCorrection, xray::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat)
    @assert isequal(unk.subshell, inner(xray)) "Unknown and X-ray don't match in XPP"
    @assert isequal(std.subshell, inner(xray)) "Standard and X-ray don't match in XPP"
    return Fχ(unk, xray, θunk) / Fχ(std, xray, θstd)
end

"""
    Z(unk::MatrixCorrection, std::MatrixCorrection)
The atomic number correction factor.
"""
function Z(unk::MatrixCorrection, std::MatrixCorrection)
    @assert isequal(unk.subshell, std.subshell)
    "Unknown and standard matrix corrections don't apply to the same sub-shell."
    return F(unk) / F(std)
end

"""
    A(unk::MatrixCorrection, std::MatrixCorrection, xray::CharXRay, χcunk=0.0, tcunk=0.0, χcstd=0.0, tcstd=0.0)
The absorption correction factors.
"""
function A(unk::MatrixCorrection, std::MatrixCorrection, xray::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat)
    @assert isequal(unk.subshell, inner(xray)) "Unknown and X-ray don't match in XPP"
    @assert isequal(std.subshell, inner(xray)) "Standard and X-ray don't match in XPP"
    return ZA(unk, std, xray, θunk, θstd) / Z(unk, std)
end



"""
    kcoating(ty::Type{<:MatrixCorrection}, subtrate::Material, coating::Material, cxr::CharXRay, e0::Real, toa::Real, τ::Real)

Estimate the k-ratio for a coating of mass-thickness τ (g/cm²) on the specified substrate. The standard for the coating
is assumed to be of the same material as the coating.
"""
function kcoating(
    ty::Type{<:MatrixCorrection},
    subtrate::Material,
    coating::Material,
    cxr::CharXRay,
    e0::Real,
    toa::Real,
    τ::Real,
)
    @assert element(cxr) in keys(coating) "$cxr must be produced by one of the elements in $(name(coating))"
    @assert e0 > energy(inner(cxr)) "The beam energy must exceed the edge energy for $cxr."
    submc, coatmc = matrixcorrection(ty, subtrate, inner(cxr), e0), matrixcorrection(ty, coating, inner(cxr), e0)
    return Fχp(submc, cxr, toa, τ) / Fχ(coatmc, cxr, toa)
end

""""
    massthickness(ty::Type{<:MatrixCorrection}, subtrate::Material, coating::Material, cxr::CharXRay, e0::Real, toa::Real, k::Real)

Estimate the mass-thickness of a ultra-thin layer of a `coating` material on a `substrate` from a measured k-ratio `k`
of a characteristic X-ray `cxr`.  Works for k-ratios of the order of 1 %.  The standard for the coating is assumed
to be of the same material as the coating.
"""
function massthickness(
    ty::Type{<:MatrixCorrection},
    substrate::Material,
    coating::Material,
    cxr::CharXRay,
    e0::Real,
    toa::Real,
    k::Real,
)
    submc, coatingmc = matrixcorrection(ty, substrate, inner(cxr), e0), matrixcorrection(ty, coating, inner(cxr), e0)
    f(τ) = k - Fχp(submc, cxr, toa, τ) / Fχ(coatingmc, cxr, toa)
    return find_zero(f, 0.01 * range(ty, substrate, e0), Roots.Order1()) / coating[element(cxr)]
end


"""
    correctkratios(krs::AbstractVector{KRatio}, coating::Material, θtoa::Real, ρz::Real)::Vector{KRatio}

This function is mainly for pedagogical purposes.  It takes a `KRatio[]`, a coating `Material` on the unknown,
and a mass-thickness (g/cm²) and creates a new `KRatio[]` that accounts for the intensity missing due to absorption
by the coating.  Favor the coating correction built into `ZAFCorrection` or `MultiZAF`.
"""
function correctkratios(krs::AbstractVector{KRatio}, coating::Material, ρz::Real)::Vector{KRatio}
    res = KRatio[]
    for kr in krs
        if !(kr.element in keys(coating))
            kratio = kr.kratio / exp(-χ(coating, brightest(kr.lines), kr.unkProps[:TakeOffAngle]) * ρz)
            push!(res, KRatio(kr.lines, kr.unkProps, kr.stdProps, kr.standard, kratio))
        else
            push!(res, kr)
        end
    end
    return res
end
