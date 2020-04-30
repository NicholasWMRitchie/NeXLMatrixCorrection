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

# Abstract class for a base fluorescence correction
include("fluorescencecorrection.jl")
export FluorescenceCorrection # Abstract class for the algorithm implementing the F-term
export NullFluorescence
export F # Compute the flurorescence correction

# Abstract class for a base coating correction algorithm
include("coating.jl")
export CoatingCorrection
export NullCoating  # 100% transmission (no correction)
export Coating      # A basic multi-layer coating for ultra-thin coatings
export transmission # Coating transmission
export carboncoating # Build a carbon coating
export coatingcorrection # Builds a coating correction for the specified algorithm (Null or Coating)

# Pulls together the MatrixCorrection, FluorescenceCorrection and CoatingCorrection
include("zafcorrection.jl")
export ZAFCorrection
export ZAF # Build a full ZAFCorrection based on your choice of algorithms
export ZAFc # Combined correction factor ZAF + coating
export matrixcorrection # Builds a
export kcoating # Calculate the k for an ultra-thin coating on a substrate
export massthickness # Estimate the mass-thickness of a coating on a substrate

# Wraps multiple ZAFCorrection objects to calculate the full matrix correction for one or more characteristic x-rays
# This addresses the need to
include("multizaf.jl")
export MultiZAF # Represents a multiline ZAF correction
export gZAFc # Combined correction factor (generation + ZAF + coating)
export k # Computed k-ratio

# Implements Pouchou & Pichoir's XPP  ϕ(ρz) model
include("xpp.jl")
export XPP # XPP structure
export xpp # Build an XPP correction

# Implements Pouchou & Pichoir's XPP  ϕ(ρz) model with full uncertainty calculation
include("xppu.jl")
export JzLabel
export Ju
export StepMJZbarb, StepDPT, StepQlaOoS

# Implements Reed's 1991 fluorescence correction
include("reed.jl")
export ReedFluorescence
export reedFluorescence # Construct a structure encapsulating the Reed fluorescence correction model

# For picking which characteristic lines to use to perform the matrix correction
include("kratioopt.jl")
export KRatioOptimizer # Abstract class
export SimpleKRatioOptimizer # A very simple implmentation of KRatioOptimizer
export optimizeks

# Performs iteration to estimate the composition from measured k-ratios
include("iterate.jl")
export UnmeasuredElementRule # Calculation elements by difference or stoichiometry or ???
export NullUnmeasuredRule # Don't do anything..
export UpdateRule # A rule implementing update(...)
export NaiveUpdateRule # Simple iteration
export WegsteinUpdateRule # Wegstein iteration (gradient)
export RecordingUpdateRule
export ConvergenceTest # Test for convergence using converged(...)
export RMSBelowTolerance, AllBelowTolerance, IsApproximate # Difference implementations of ConvergenceTest
export Iteration # Defines the iteration procedure
export IterationResult # The output from iterateks(...)
export iterateks # Perform the iteration
export quantify # Exported when NeXLSpectrum is loaded

function __init__()
    @require Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004" include("gadflysupport.jl")
    @require NeXLSpectrum = "6c578565-ca7f-4012-afc4-b8412d85af92" include("spectrumsupport.jl")
end

end
