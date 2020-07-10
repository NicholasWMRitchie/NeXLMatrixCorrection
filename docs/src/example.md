## A Simple Example using NeXLMatrixCorrection
Load the necessary libraries
````julia
using NeXLMatrixCorrection  # Inplements the algorithms (auto-loads NeXLCore)
using DataFrames            # For the table
````




##### Convert k-ratios to composition.
Define the measurement conditions
````julia
lbl = label("K458")  # This labels the measurement
# Define the measurement properties (:BeamEnery and :TakeOffAngle are required by most matrix correction algorithms)
unkProps = Dict(:BeamEnergy=>15.0e3, :TakeOffAngle=>deg2rad(40.0), :Coating=>Film(pure(n"C"), 7.0e-7))
stdProps = Dict(:BeamEnergy=>15.0e3, :TakeOffAngle=>deg2rad(40.0), :Coating=>Film(pure(n"C"), 15.0e-7))
# Create a list of the measurement k-ratios.
krs = [
    KRatio([n"O K-L3"], unkProps, stdProps, mat"SiO2", uv(0.746227,0.0010) ),
    KRatio([n"Si K-L3"], unkProps, stdProps, mat"SiO2", uv(0.441263,0.0012) ),
    KRatio([n"Zn K-L3"], unkProps, stdProps, mat"Zn", uv(0.027776,0.0002) ),
    KRatio([n"Ba L3-M5"], unkProps, stdProps, mat"BaCl", uv(0.447794,0.0020) )
]
````


````
4-element Array{KRatio,1}:
 k[SiO2, O K-L3] = 0.74623 ± 0.001
 k[SiO2, Si K-L3] = 0.44126 ± 0.0012
 k[Zn, Zn K-L3] = 0.027776 ± 0.0002
 k[BaCl, Ba L3-M5] = 0.44779 ± 0.002
````




##### Perform the Iteration
````julia
# Now perform the iteration on the k-ratios
res = quantify(lbl, krs)
# Tabulate the results...
asa(DataFrame, res, withZAF=true)
````


````
4×14 DataFrame. Omitted printing of 8 columns
│ Row │ Label │ Element │ Standard │ Lines    │ Mass Frac. │ Δ[Mass Frac.] 
│
│     │ Label │ String  │ String   │ String   │ Float64    │ Float64       
│
├─────┼───────┼─────────┼──────────┼──────────┼────────────┼───────────────
┤
│ 1   │ K458  │ Si      │ SiO2     │ Si K-L3  │ 0.23133    │ 0.00063       
│
│ 2   │ K458  │ Ba      │ BaCl     │ Ba L3-M5 │ 0.4218     │ 0.0019        
│
│ 3   │ K458  │ O       │ SiO2     │ O K-L3   │ 0.32701    │ 0.00044       
│
│ 4   │ K458  │ Zn      │ Zn       │ Zn K-L3  │ 0.030727   │ 0.00022       
│
````





Now take a slightly different perspective that focuses more on iteration related data.
````julia
asa(DataFrame, res, withZAF=false)
````


````
4×8 DataFrame. Omitted printing of 2 columns
│ Row │ Label │ Element │ Converged │ Iterations │ Mass Frac. │ Δ[Mass Frac
.] │
│     │ Label │ String  │ Bool      │ Int64      │ Float64    │ Float64    
   │
├─────┼───────┼─────────┼───────────┼────────────┼────────────┼────────────
───┤
│ 1   │ K458  │ Si      │ 1         │ 7          │ 0.23133    │ 0.00063    
   │
│ 2   │ K458  │ Ba      │ 1         │ 7          │ 0.4218     │ 0.0019     
   │
│ 3   │ K458  │ O       │ 1         │ 7          │ 0.32701    │ 0.00044    
   │
│ 4   │ K458  │ Zn      │ 1         │ 7          │ 0.030727   │ 0.00022    
   │
````





Nicholas W. M. Ritchie, 30-Apr-2020
