"""
    NullCorrection
Implements Castaing's First Approximation (i.e. No correction Z⋅A = 1)
"""
struct NullCorrection <: MatrixCorrection
    material::Material
    subshell::AtomicSubShell
    E0::Float64

    NullCorrection(mat::Material, ashell::AtomicSubShell, e0) = new(mat, ashell, e0)
end

Base.show(io::IO, nc::NullCorrection) = print(io, "Null[ $(nc.material), $(nc.subshell), $(0.001 * nc.E0) keV]")

ℱ(mc::NullCorrection) = 1.0
ℱχ(mc::NullCorrection, χ::AbstractFloat) = 1.0
NeXLCore.atomicsubshell(mc::NullCorrection) = mc.subshell
NeXLCore.material(mc::NullCorrection) = mc.material
beamEnergy(mc::NullCorrection) = mc.E0 # in eV
ϕ(::NullCorrection, ρz) = @assert "Does not make sense for the null correction."
ϕabs(::MatrixCorrection, ρz, χ) = @assert "Does not make sense for the null correction."
NeXLMatrixCorrection.matrixcorrection(::Type{NullCorrection}, mat::Material, ass::AtomicSubShell, e0::AbstractFloat) = NullCorrection(mat,ass,e0)
NeXLCore.minproperties(::Type{NullCorrection}) = ( )
