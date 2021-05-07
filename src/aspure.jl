
"""
aspure(krs::Union{KRatio, KRatios}, mc::Type{<:MatrixCorrection} = XPP, fc::Type{<:FluorescenceCorrection} = NullFluorescence, cc::Type{<:CoatingCorrection} = NullCoating)

Reinterprete a k-ratio or k-ratios as relative to a pure, uncoated element at the same beam energy and take-off angle as the unknown.
"""
function aspure(
    krs::KRatio, #
    mc::Type{<:MatrixCorrection} = XPP,
    fc::Type{<:FluorescenceCorrection} = NullFluorescence,
    cc::Type{<:CoatingCorrection} = NullCoating,
)
    res = krs
    if (length(keys(krs.standard)) > 1) || (krs.unkProps[:TakeOffAngle]!=krs.stdProps[:TakeOffAngle]) || (krs.unkProps[:BeamEnergy]!=krs.stdProps[:BeamEnergy])
        pureprops = Dict( :TakeOffAngle => krs.unkProps[:TakeOffAngle], :BeamEnergy => krs.unkProps[:BeamEnergy] )
        purestd = NeXLCore.pure(krs.element)
        pure = zafcorrection(mc, fc, cc, purestd, krs.xrays, pureprops[:BeamEnergy], missing)
        std = zafcorrection(mc, fc, cc, krs.standard, krs.xrays, krs.stdProps[:BeamEnergy], get(krs.stdProps, :Coating, missing))        
        kpure = k(pure, std, pureprops[:TakeOffAngle], krs.stdProps[:TakeOffAngle])
        res = KRatio(krs.xrays, krs.unkProps, pureprops, purestd, krs.kratio / kpure)
    end
    return res
end
function aspure(
    krs::KRatios, #
    mc::Type{<:MatrixCorrection} = XPP,
    fc::Type{<:FluorescenceCorrection} = NullFluorescence,
    cc::Type{<:CoatingCorrection} = NullCoating,
)
    res = krs
    if (length(keys(krs.standard)) > 1) || (krs.unkProps[:TakeOffAngle]!=krs.stdProps[:TakeOffAngle]) || (krs.unkProps[:BeamEnergy]!=krs.stdProps[:BeamEnergy])
        pureprops = Dict( :TakeOffAngle => krs.unkProps[:TakeOffAngle], :BeamEnergy => krs.unkProps[:BeamEnergy] )
        purestd = NeXLCore.pure(krs.element)
        pure = zafcorrection(mc, fc, cc, purestd, krs.xrays, pureprops[:BeamEnergy], missing)
        std = zafcorrection(mc, fc, cc, krs.standard, krs.xrays, krs.stdProps[:BeamEnergy], get(krs.stdProps, :Coating, missing))        
        kpure = k(pure, std, pureprops[:TakeOffAngle], krs.stdProps[:TakeOffAngle])
        res = KRatios(krs.xrays, krs.unkProps, pureprops, purestd, krs.kratios / kpure)
    end
    return res
end