"""
    thinfilm(krs::Vector{KRatio})::Tuple{Float64,Composition}

Performs a single-layer thin-film correction on the supplied k-ratio data.
Returns a tuple containing the estimated thickness and composition.
"""
function thinfilm(krs::Vector{KRatio})::Tuple{Float64,Composition}
   function ZAFc(std::MultiZAF, θstd::AbstractFloat)
      a, n = 0.0, 0.0
      for (sh, cxrs2) in splitbyshell(std)
         zafS = std.zafs[sh]
         for cxr in cxrs2
            w = weight(NormalizeByShell, cxr)
            a += w * ℱχ(zafS.za, cxr, θstd) * F(zafS.f, cxr, θstd) * transmission(zafS.coating, cxr, θstd)
            n += w
         end
      end
      return a / n
   end
   function absorb(comp, cxrs, ρz, toa)
      if ρz > 0
         a, n = 0.0, 0.0
         for cxr in cxrs
            w = weight(NormalizeByShell, cxr)
            χ = mac(comp, cxr) / sin(toa)
            a += w * (1.0 - Math.exp(-χ * ρz)) / (χ * ρz)
            n += w
         end
         a/n
      else 
         1.0
      end
   end
   ρz, comp = -1000.0, nothing
   for _ in 1:10
      ρz2 = 0.0
      mfs = HashMap{Element,Float64}()
      for kr in krs
         toa, e0 = kr.stdProps[:TakeOffAngle], kr.stdProps[:BeamEnergy]
         std, xrs, el = standard(kr), xrays(kr), element(kr)
         xpp = zafcorrection(XPP, ReedFluorescence, NullCoating, std, xrs, e0, nothing)
         # f is the absorption correction factor
         tmp = (std[el] * ZAFc(xpp, toa) / absorb(std, xrs, ρz, toa)) * kratio(kr)
         mfsρz[el] = tmp
         ρz2 += tmp
      end
      mfs = map((k, v) -> k => v / ρz2, mfsρz)
      comp2 = Material("thin-film", mfs)
      breakout = Math.abs(rhoz - rhoz2) / rhoz2 < 0.005
      comp, ρz = comp2, rhoz2
      if breakout
         break
      end
   end
   return tuple(ρz, comp)
end


function multilayer(krs::Vector{KRatio}, layer::AbstractMap{Element,Integer})::Vector{Pair{Material, Float64}}
   return Pair{Material, Float64}[]
end
