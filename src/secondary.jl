
"""
   lenardCoefficient(e0::Float64, ashell::AtomicShell)

Computes the Lenard coefficient according to the algorithm of Heinrich.
"Heinrich K. F. J. (1967) EPASA 2, paper no. 7"
"""
lenardCoefficient(e0::Float64, ashell::AtomicShell) =
   4.5e5 / ((0.001*e0)^1.65 - (0.001*energy(ashell))^1.65)

"""
    ionizationDepthRatio(primary::AtomicShell, secondary::AtomicShell, e0::Float64)

Ionization depth ratio from "Reed S.J.B. (1990) Microbeam Analysis, p.109"
"""
function ionizationDepthRatio(primary::AtomicShell, secondary::AtomicShell, e0::Float64)
   uA, uB = e0 / energy(secondary), e0 / energy(primary)
   return ((((uB * log(uB)) - uB) + 1.0) / (((uA * log(uA)) - uA) + 1.0))
end

"""
    familyFactor(shellA::AtomicShell, shellB::AtomicShell)::Float64
Accounts for the differences in ionization cross section between K , L and M shells
"""
function familyFactor(shellA::AtomicShell, shellB::AtomicShell)::Float64
   fA, fB = family(shellA), family(shellB)
   if fA == fB
      res = 1.0 # fA==fB
   else
      if (fA == 'K') && (fB == 'L')
         res = 1.0 / 0.24
      elseif (fA == 'L') && (fB == 'K')
         res = 0.24
      elseif fA == 'M'
         @assert((fB == 'K') || (fB == 'L'))
         res = 0.02
      else
         res = 0.0
      end
   end
   return res
end

"""
    ionizationFraction(shell::AtomicShell)

The fraction of the ionizations to attribute to the specified shell.  Computed
from the jump ratio.
"""
function ionizationFraction(shell::AtomicShell)
   r = jumpRatio(shell)
   return r >= 1.0 ? (r - 1.0) / r : 0.0
end

struct ReedInternal
   primary::CharXRay
   kk::Float64
   lenard::Float64
   muB::Float64

   function ReedInternal(comp::Material, primary::CharXRay, secondary::AtomicShell, e0::Float64)
      bElm, aElm = element(primary), element(secondary)
      @assert (aElm in keys(comp)) && (bElm in keys(comp))
      cB = comp[bElm]
      muB_A, muB = mac(aElm, primary), mac(comp, primary)
      ionizeF = ionizationFraction(secondary)
      fluorB = meanFluorescenceYield(inner(primary))
      v = lenardCoefficient(e0, secondary) / muB # keV
      ss = ionizationDepthRatio(inner(primary), secondary, e0)
      f = familyFactor(secondary, inner(primary));
      k = f * 0.5 * cB * (muB_A / muB) * ionizeF * fluorB * (a(aElm) / a(bElm)) * ss
      return new(primary, k, v, muB)
   end
end

Base.show(io::IO, ri::ReedInternal) =
   print(io,repr(ri.primary))

function Finternal(ri::ReedInternal, secondary::CharXRay, toa::Float64, comp::Material)
   @assert(ri.kk > 0.0)
   u = mac(comp, secondary) / (sin(toa) * ri.muB)
   # TODO: Evaluate whether weight(ri.primary) is necessary/correct???
   return weight(ri.primary) * ri.kk * ((log(1.0 + u) / u) + (log(1.0 + ri.lenard) / ri.lenard))
end

struct ReedFluorescence <: FluorescenceCorrection
   comp::Material
   secondary::AtomicShell # Emitter
   e0::Float64 # Incident beam energy
   # Cached partial calculations
   exciters::Vector{ReedInternal}
end

Base.show(io::IO, reed::ReedFluorescence) =
   print(io,"Reed[$(reed.secondary) due to $(reed.exciters) for $(name(reed.comp)) at $(reed.e0/1000.0) keV]")

"""
   F(reed::ReedFluorescence, secondary::CharXRay, toa::Float64)

Compute the enhancement of the secondary characteristic X-ray due to the
primaries specified in reed.
"""
F(reed::ReedFluorescence, secondary::CharXRay, toa::Float64) =
   mapreduce(ex -> Finternal(ex, secondary, toa, reed.comp), +, reed.exciters)

"""
    reedFluorescence(comp::Material, primary::Vector{CharXRay}, secondary::AtomicShell, e0::Float64)

Construct an instance of a ReedFluorescence correction structure to compute the
secondary fluorescence due to a primary characteristic X-ray in the specified
material and beam energy.
"""
function reedFluorescence(
   comp::Material,
   primarys::Vector{CharXRay},
   secondary::AtomicShell,
   e0::Float64,
)
   ris = Vector{ReedInternal}()
   if element(secondary) in keys(comp)
      for primary in primarys
         if (energy(primary) >= energy(secondary)) &&
            (element(primary) in keys(comp))
            push!(ris, ReedInternal(comp, primary, secondary, e0))
         end
      end
   end
   return ReedFluorescence(comp, secondary, e0, ris)
end

"""
    reedFluorescence(comp::Material, secondary::AtomicShell, e0::Float64)

Construct an instance of a ReedFluorescence correction structure to compute the
secondary fluorescence in the specified material and beam energy.
"""
function reedFluorescence(
   comp::Material,
   secondary::AtomicShell,
   e0::Float64;
   eThresh = 5.0e3,
   wThresh = 0.01,
)
   inERange(cxr::CharXRay, ee::Float64) =
      (energy(cxr) > ee) && (energy(cxr) < ee + eThresh)
   primaries = Vector{CharXRay}()
   for elm in keys(comp)
      candidates = characteristic(elm, alltransitions, wThresh / comp[elm], e0)
      for cxr in filter(c -> inERange(c, energy(secondary)), candidates)
         push!(primaries, cxr)
      end
   end
   return reedFluorescence(comp, primaries, secondary, e0)
end
