
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

struct ReedFluorescence <: FluorescenceCorrection
   comp::Material
   primary::CharXRay # Exciter
   secondary::AtomicShell # Emitter
   e0::Float64 # Incident beam energy

   kk::Float64  # Overall constant factor
   lenard::Float64 # Lenard coefficient
   muB::Float64 # MAC for primary in comp
end

Base.show(io::IO, reed::ReedFluorescence) =
   print(io,"Reed[$(reed.secondary) due to $(reed.primary) for $(name(reed.comp)) at $(reed.e0/1000.0) keV]")

"""
   F(reed::ReedFluorescence, secondary::CharXRay, toa::Float64)

Compute the enhancement of the secondary characteristic X-ray due to the
primary specified in reed.
"""
function F(reed::ReedFluorescence, secondary::CharXRay, toa::Float64)
   if reed.kk > 0.0
      u = mac(reed.comp, secondary) / (sin(toa) * reed.muB)
      return reed.kk * ((log(1.0 + u) / u) + (log(1.0 + reed.lenard) / reed.lenard))
   end
   return 0.0
end

"""
    reedFluorescence(comp::Material, primary::CharXRay, secondary::AtomicShell, e0::Float64)

Construct an instance of a ReedFluorescence correction structure to compute the
secondary fluorescence due to a primary characteristic X-ray in the specified
material and beam energy.
"""
function reedFluorescence(comp::Material, primary::CharXRay, secondary::AtomicShell, e0::Float64)
   if energy(primary) >= energy(secondary)
      bElm, aElm = element(primary), element(secondary)
      if (aElm in keys(comp)) && (bElm in keys(comp))
         cB = comp[bElm]
         muB_A, muB = mac(aElm, primary), mac(comp, primary)
         ionizeF = ionizationFraction(secondary)
         fluorB = meanFluorescenceYield(inner(primary))
         v = lenardCoefficient(e0, secondary) / muB # keV
         ss = ionizationDepthRatio(inner(primary), secondary, e0)
         f = familyFactor(secondary, inner(primary));
         k = f * 0.5 * cB * (muB_A / muB) * ionizeF * fluorB * (a(aElm) / a(bElm)) * ss
         return ReedFluorescence(comp, primary, secondary, e0, k, v, muB)
      end
      return ReedFluorescence(comp, primary, secondary, e0, 0.0, 0.0, 0.0)
   end
end
