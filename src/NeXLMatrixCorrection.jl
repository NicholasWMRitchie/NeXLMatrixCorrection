module NeXLMatrixCorrection

using Reexport
using Requires

@reexport using NeXLCore  # Re-export so NeXLCore automatically is available with being explicitly loaded

# Abstract class for a base matrix correction algorithm (like XPP or ???)
include("matrixcorrection.jl")
export MatrixCorrection # Abstract class for the algorithm implementing the Z- & A-terms
export Z # Atomic number correction
export A # Absorption correction
export ZA # Combined ZA (ϕ(ρz)) correction (more efficient than Z(...)*A(...))
export ϕ # ϕ(ρz) function
export ϕabs # ϕ(ρz) function (absorbed)
export F  # Generated intensity function
export Fχ # Absorbed intensity function
export matrixcorrection # Factory method for MatrixCorrection algorithms (XPP and NullCorrection)
export continuumcorrection # Factory method for creating matrix correction algorithms for continuum correction
export correctcontinuum # The method to calculate the continuum correction
export χ # Calculates reduced mass-absorption coefficient

# Abstract class for a base fluorescence correction
include("fluorescencecorrection.jl")
export FluorescenceCorrection # Abstract class for the algorithm implementing the F-term
export NullFluorescence
# export F # Compute the flurorescence correction
export fluorescencecorrection # Factory method for FluorescenceCorrection algorithms (Reed and NullFluorescence)

# Abstract class for a base coating correction algorithm
include("coating.jl")
export CoatingCorrection
export NullCoating  # 100% transmission (no correction)
export Coating      # A basic multi-layer coating for ultra-thin coatings
export carboncoating # Build a carbon coating
export coatingcorrection # Factory method for CoatingCorrection algorithms (Null or Coating)

# Pulls together the MatrixCorrection, FluorescenceCorrection and CoatingCorrection into a single object
include("zafcorrection.jl")
export ZAFCorrection
export zafcorrection # Factory method for ZAFCorrection
export ZAFc # Combined correction factor ZAF + coating
export beamEnergy

# Use massthickness to estimate the thickness of an ultra-thin coating on a substrate.
export massthickness # Estimate the mass-thickness of a coating on a substrate
export coatingasfilm # Calculates the coating thickness and creates a Film object to represent it.
export kcoating # Calculate the measured k for an ultra-thin coating on a substrate

# Wraps multiple ZAFCorrection objects to calculate the full matrix correction for one or more characteristic x-rays
# This addresses the need to
include("multizaf.jl")
export MultiZAF # Represents a multiline ZAF correction
export gZAFc # Combined correction factor (generation + ZAF + coating)
export k # Computed k-ratio
export aspure # The k-ratio for this X-ray relative to a pure element (use with NeXLSpectrum.VectorQuant)

# Implements Pouchou & Pichoir's XPP  ϕ(ρz) model
include("xpp.jl")
export XPP # <: MatrixCorrection
export xpp # Build an XPP correction

include("citzaf.jl")
export CitZAF # <: MatrixCorrection

include("riveros.jl")
export Riveros1993 # <: MatrixCorrection

# Implements Pouchou & Pichoir's XPP  ϕ(ρz) model with full uncertainty calculation
include("xppu.jl")
export JzLabel
export Ju
export StepMJZbarb, StepDPT, StepQlaOoS

# Implements Merlet's XPhi algorithm
include("xphi.jl")
export XPhi

# Implements the null matrix correction
include("nullmatrixcorrection.jl")
export NullCorrection

# Implements Reed's 1991 fluorescence correction
include("reed.jl")
export ReedFluorescence
export reedFluorescence # Construct a structure encapsulating the Reed fluorescence correction model

# For picking which characteristic lines to use to perform the matrix correction
include("kratioopt.jl")
export KRatioOptimizer # Abstract class
export SimpleKRatioOptimizer # A very simple implmentation of KRatioOptimizer
export optimizeks # The method required of KRatioOptimizer

include("standard.jl")
export Standard # Collects k-ratios to be used for standardization.
export standardize # Apply similar standards to a KRatio

# Performs iteration to estimate the composition from measured k-ratios
include("iterate.jl")
# Abstract class for calculating elements by difference or stoichiometry or ???
export UnmeasuredElementRule
export NullUnmeasuredRule # Don't do anything..

# Abstract class for updating between iteration steps
export UpdateRule # A rule implementing update(...)
export NaiveUpdateRule # Simple iteration
export WegsteinUpdateRule # Wegstein iteration (gradient)
export RecordingUpdateRule

export ConvergenceTest # Test for convergence using converged(...)
export RMSBelowTolerance, AllBelowTolerance, IsApproximate # Difference implementations of ConvergenceTest

export Iteration # Defines the iteration procedure
export IterationResult # The output from quantify(...)

export quantify # Perform the iteration on KRatio(s) or FilterFitResult

include("supportedthinfilms.jl")
export SupportedThinFilms
export nlayers
export χs
export outer, inner

#include("xfilm.jl")
#export XFilm

include("helpers.jl")
export zaf # Generates a table of ZAF correction factors

include("defaultstandards.jl")
export getstandards # Get a list of suggested standards for an element.

function __init__()
    @require Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004" include("gadflysupport.jl")
end

end
