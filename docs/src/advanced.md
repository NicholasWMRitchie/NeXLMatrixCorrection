# Advanced Matrix Correction

The matrix correction algorithm is designed to permit certain advanced options 
as alternative mechanisms for computing "unmeasured elements", alternative 
iteration algorithms and alternative halt algorithms.

## Unmeasured Elements
It is not always possible or desirable to measure all the elements present 
in a sample.  Regardless, all matrix correction algorithms need to know the 
mass fraction of every element in the material since all the elements 
contribute to the stopping power, backscatter correction, absorption
and other corrections.

One way to get around this requirement is to provide alternative mechanisms 
to compute the unmeasured element or elements.  The classic example of an 
unmeasured element algorithm is "oxygen by stoichiometry."  This algorithm 
assumes that each non-oxygen element is associated with a fixed number of 
oxygen atoms.  Thus commonly it is assumed that each Si atom is associated 
with two O atoms and two Al atoms are associated with three O atoms.

Oxygen-by-stoichiometry is implemented in `NeXLMatrixCorrection` by 
implementing the `UnmeasuredElementRule` abstract type. Types implementing 
this type are expected to implement the methods:

```julia
    NeXLUncertainties.compute(
        uer::<:UnmeasuredElementRule, 
        inp::Dict{Element,<:AbstractFloat}
    )::Dict{Element,AbstractFloat}
    
    isunmeasured(uer::<:UnmeasuredElementRule, elm::Element)::Bool
```

The first method is very general. Given one input mapping from element to mass-fraction,
it computes a second output mapping from element to mass-fraction.  Typically, the input 
mass-fractions are replicated in the output and additional elements and associated 
mass-fractions are added.  The second method serves to identify which elements are
computed by the `UnmeasuredElementRule`.

`NeXLMatrixCorrection` implements `ElementByStoichiometry` which is a generalization
of oxygen-by-stoichiometry and `ElementByFiat` which simply adds a fixed mass-fraction
of a specified element.

The method `OByStoichiometry(valences = NeXLCore.defaultValences)` creates an 
`ElementByStoichiometry` with the specified valences.  To specify custom valences,
copy the default set and update the copy with your desired values. 

To sequentially apply a set of `UnmeasuredElementRule`s use the 
`MultiUnmeasuredElementRule` type.

You can implement custom `UnmeasuredElementRule`s to perform other modes.

## Iteration algorithms

Iteration is implemented by concrete types implementing the `UpdateRule`
abstract type.   By default, `NeXLMatrixCorrection` uses the Wegstein iteration 
algorithm as implemented by `WegsteinUpdateRule`.  This is an efficient and
effective iteration algorith.  Alternatively, the `NaiveUpdateRule` is also 
available.  Other rules can be implemented by adding a concreate type extending
`UpdateRule` and implementing the method:

```julia
update(
    ::Type{<:UpdateRule},
    prevcomp::Material,
    measured::Vector{KRatio},
    zafs::Dict{Element, AbstractFloat}
)::Dict{Element, Float64}  
```

This rule takes a previous estimate of the composition, a measured vector of `KRatio`s and
a mapping of matrix correction factors to estimate the mass fraction of all **measured** 
elements.  Any `UnmeasuredElementRule`s are applied subsequently to fully estimate the
matrix composition.

### Convergence Tests

An additional abstract type `ConvergenceTest` is provided for alternative mechanisms to indicate
that the iteration should be terminated.

Data types based on `ConvergenceTest` should implement the function
```julia
converged(
    ia::ConvergenceTest, 
    meas::Vector{KRatio}, 
    computed::Dict{Element,<:AbstractFloat}
)::Bool
```
This function returns true when the computed k-ratios match the measured k-ratios to within a 
tolerance.

### Configuring Iteration

The `quantify(...)` function is configured with an `Iteration` data item.  The `Iteration` item 
specifies the `MatrixCorrection`, `FluorescenceCorrection` and `CoatingCorrection` algorithms as 
well as the `UpdateRule`, the `ConvergenceTest` and the `UnmeasuredElementRule`.
```julia
Iteration(
        mct::Type{<:MatrixCorrection},
        fct::Type{<:FluorescenceCorrection},
        cct::Type{<:CoatingCorrection};
        updater = WegsteinUpdateRule(),
        converged = RMSBelowTolerance(0.00001),
        unmeasured = NullUnmeasuredRule(),
    )
```
The `quantify(...)` function takes the `Iteration` item and the initial data.  The `coating`
argument is used to compute the coating thickness from the measured k-ratios.  It takes a 
`Pair{CharXRay, Material}` which defines the X-ray to use to estimate the coating thickness and
the coating material.  The X-ray must only be produced in the coating and not the sample for this 
mechanism to work.

```julia
quantify(
    name::Union{String, Label},
    measured::Vector{KRatio},
    iteration::Iteration = Iteration(XPP, ReedFluorescence, Coating);
    maxIter::Int = 100, 
    estComp::Union{Nothing,Material}=nothing, 
    coating::Union{Nothing, Pair{CharXRay, <:Material}}=nothing
)::IterationResult

quantify(
    measured::AbstractVector{KRatios{T}},
    iteration::Iteration = Iteration(XPP, NullFluorescence, NoCoating);
    kro::KRatioOptimizer = SimpleKRatioOptimizer(1.5),
    maxIter::Int = 100, 
    estComp::Union{Nothing,Material}=nothing, 
    coating::Union{Nothing, Pair{CharXRay, <:Material}}=nothing
)
```

There are also various specializations to handle fitted spectrum data and
spectrum images.  These are listed in the online documentation for various 
other packages.