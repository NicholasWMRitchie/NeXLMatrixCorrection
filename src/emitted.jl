"""
    emitted_intensities(comp::Material, elm::Element, dose::Float64, e0::Float64, θtoa::Float64, Ω::Float64; mc=XPP, fc=ReedFluorescence)
    emitted_intensities(comp::Material, dose::Float64, e0::Float64, θtoa::Float64, Ω::Float64; mc=XPP, fc=ReedFluorescence)

  * comp: The composition of the Material
  * elm: The ionized element
  * dose: Electron dose in nA⋅s
  * e0: Beam energy in eV
  * θtoa: Take off angle in radians
  * Ω: Solid angle steradians
  
  * returns Dict{CharXRay, Float64} containing characteristic X-rays and the emitted intensities.

Computes the intensity emitted from the sample at the take-off angle into the specified solid angle.
The XPP and Reed fluorescence algorithms are used by default but any ϕ(ρz)-model in which the
integral Fχ is normalized to the emission from a single shell may be used.
"""
function emitted_intensities(comp::Material, elm::Element, dose::Float64, e0::Float64, θtoa::Float64, Ω::Float64; mc=XPP, fc=ReedFluorescence, minweight=0.001)
    cxrs = characteristic(elm, alltransitions, minweight, e0)
    kk = atoms_per_g(comp, elm) *     # Number of atoms of elm in 1 g of comp
        6.241509074460764e9 * dose *  # ustrip(NoUnits, (dose*u"nA*s")/ElementaryCharge) * # Number of electrons per nA⋅s
        Ω/(4π)                        # solid angle
    si = Dict{CharXRay, Float64}()
    for ash in unique(inner.(cxrs))
        zaf = zafcorrection(mc, fc, NullCoating, comp, ash, e0, missing)
        icx = ionizationcrosssection(ash, e0) * occupancy(ash) # Cross section for ionization of ash per atom 
        for cxr in cxrs
            ashi, inn, out = NeXLCore.subshellindex(ash), NeXLCore.innerindex(cxr), NeXLCore.outerindex(cxr)
            if ashi <= inn
                fy = NeXLCore.xrayweight(NormalizeRaw, cxr.z, ashi, inn, out)
                if fy > 0.0
                    si[cxr] = get(si, cxr, 0.0) +
                        kk * icx *              # constants 
                        fy *                    # Fraction of ionizations that eventually relax via cxr
                        Fχ(zaf.za, cxr, θtoa) * # Absorption corrected integral of ϕ(ρz) 
                        F(zaf.f, cxr, θtoa)     # Fluorescence correction factor (~1.)
                end
            end
        end
    end
    return si
end
function emitted_intensities(comp::Material, dose::Float64, e0::Float64, θtoa::Float64, Ω::Float64; kwargs...)
    mapreduce(merge, elms(comp)) do elm
        emitted_intensities(comp, elm, dose, e0, θtoa, Ω; kwargs...)
    end
end
