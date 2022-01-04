
"""
   lenardcoefficient(e0::Float64, ashell::AtomicSubShell)

Computes the Lenard coefficient according to the algorithm of Heinrich.
"Heinrich K. F. J. (1967) EPASA 2, paper no. 7"
"""
lenardcoefficient(e0::Float64, ashell::AtomicSubShell) = 4.5e5 / ((0.001 * e0)^1.65 - (0.001 * energy(ashell))^1.65)

"""
    ionizationdepthratio(primary::AtomicSubShell, secondary::AtomicSubShell, e0::Float64)

Ionization depth ratio from "Reed S.J.B. (1990) Microbeam Analysis, p.109"
"""
function ionizationdepthratio(primary::AtomicSubShell, secondary::AtomicSubShell, e0::Float64)
    uA, uB = e0 / energy(secondary), e0 / energy(primary)
    return ((((uB * log(uB)) - uB) + 1.0) / (((uA * log(uA)) - uA) + 1.0))
end

"""
    familyfactor(shellA::AtomicSubShell, shellB::AtomicSubShell)::Float64

Accounts for the differences in ionization cross section between K , L and M shells
"""
function familyfactor(shellA::AtomicSubShell, shellB::AtomicSubShell)::Float64
    fA, fB = shell(shellA), shell(shellB)
    if fA == fB
        res = 1.0 # fA==fB
    else
        if (n(fA) == 1) && (n(fB) == 2)
            res = 1.0 / 0.24
        elseif (n(fA) == 2) && (n(fB) == 1)
            res = 0.24
        elseif (n(fA) == 3) && ((n(fB) == 1) || (n(fB) == 2))
            res = 0.02
        else
            res = 0.0
        end
    end
    return res
end

"""
    ionizationfraction(shell::AtomicSubShell)

The fraction of the ionizations to attribute to the specified shell.  Computed
from the jump ratio.
"""
function ionizationfraction(shell::AtomicSubShell)
    r = jumpratio(shell)
    return r >= 1.0 ? (r - 1.0) / r : 0.0
end

"""
The `ReedInternal` structure optimizes the calculation of the Reed correction algorithm.
"""
struct ReedInternal
    primary::CharXRay
    kk::Float64
    lenard::Float64
    muB::Float64

    function ReedInternal(comp::Material, primary::CharXRay, secondary::AtomicSubShell, e0::Float64)
        bElm, aElm = element(primary), element(secondary)
        @assert aElm in keys(comp) "Element $(aElm.symbol) not in $(comp). (aElm)"
        @assert bElm in keys(comp) "Element $(bElm.symbol) not in $(comp). (bElm)"
        cB = comp[bElm]
        muB_A, muB = mac(aElm, primary), mac(comp, primary)
        ionizeF = ionizationfraction(secondary)
        fluorB = meanfluorescenceyield(element(primary), shell(inner(primary)))
        v = lenardcoefficient(e0, secondary) / muB # keV
        ss = ionizationdepthratio(inner(primary), secondary, e0)
        f = familyfactor(secondary, inner(primary))
        k = f * 0.5 * cB * (muB_A / muB) * ionizeF * fluorB * (a(aElm) / a(bElm)) * ss
        @assert k >= 0.0 "k<0 in RI[$comp, $primary, $secondary, $e0] - $k"
        return new(primary, k, v, muB)
    end
end

"""
The `ReedFluorescence` structure implements `FluorescenceCorrection` for the Reed fluorescence model.
"""
struct ReedFluorescence <: FluorescenceCorrection
    comp::Material
    secondary::AtomicSubShell # Emitter
    e0::Float64 # Incident beam energy
    # Cached partial calculations
    exciters::Vector{ReedInternal}
end

Base.show(io::IO, reed::ReedFluorescence) =
    print(io, "Reed[$(reed.secondary) due to $(reed.exciters) for $(name(reed.comp)) at $(reed.e0/1000.0) keV]")

"""
   F(reed::ReedFluorescence, secondary::CharXRay, toa::Float64)

Compute the enhancement of the secondary characteristic X-ray due to the
primaries specified in reed.
"""
function F(reed::ReedFluorescence, secondary::CharXRay, toa::Float64)
    function finternal(ri::ReedInternal, secondary::CharXRay, toa::Float64, comp::Material)
        @assert ri.kk >= 0.0 "ri.kk = $(ri.kk) for $(ri) and $(secondary) in $(comp) - $ri"
        u = max(mac(comp, secondary) / (sin(toa) * ri.muB), 1.0e-6)
        # TODO: Evaluate whether weight(ri.primary) is necessary/correct???
        return weight(NormalizeByShell, ri.primary) * ri.kk * ((log(1.0 + u) / u) + (log(1.0 + ri.lenard) / ri.lenard))
    end
    return isempty(reed.exciters) ? 1.0 :
           1.0 + sum(ex -> finternal(ex, secondary, toa, reed.comp), reed.exciters)
end
"""
    fluorescencecorrection(::Type{ReedFluorescence}, comp::Material, primary::Vector{CharXRay}, secondary::AtomicSubShell, e0::Float64)

Construct an instance of a ReedFluorescence correction structure to compute the
secondary fluorescence due to a primary characteristic X-ray in the specified
material and beam energy.
"""
function fluorescencecorrection(
    ::Type{ReedFluorescence},
    comp::Material,
    primarys::Vector{CharXRay},
    secondary::AtomicSubShell,
    e0::Float64,
)
    ris = Vector{ReedInternal}()
    if element(secondary) in keys(comp)
        for primary in filter(p->(energy(p) >= energy(secondary)) && (element(p) in keys(comp)), primarys)
            push!(ris, ReedInternal(comp, primary, secondary, e0))
        end
    end
    return ReedFluorescence(comp, secondary, e0, ris)
end

NeXLCore.minproperties(::Type{ReedFluorescence}) = (:BeamEnergy, :TakeOffAngle)
