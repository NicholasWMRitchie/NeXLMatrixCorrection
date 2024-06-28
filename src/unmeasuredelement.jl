"""
The `UnmeasuredElementRule` mechanism provides a method to implement rules for adding unmeasured elements to
the fitting process.  Examples include element-by-stoichiometry or element-by-difference.

`UnmeasuredElementRule` implementations must include methods for:

    NeXLUncertainties.compute(uer::<:UnmeasuredElementRule, inp::Dict{Element,<:AbstractFloat})::Dict{Element,AbstractFloat}
    isunmeasured(uer::<:UnmeasuredElementRule, elm::Element)::Bool
"""
abstract type UnmeasuredElementRule end

"""
The NullUnmeasuredRule adds no additional elements in the iteration process.
"""
struct NullUnmeasuredRule <: UnmeasuredElementRule end

"""
    NeXLUncertainties.compute(::NullUnmeasuredRule, inp::Dict{Element,<:AbstractFloat})::Dict{Element,AbstractFloat}

A null UnmeasuredElementRule.  Just returns the inputs.
"""
NeXLUncertainties.compute(::NullUnmeasuredRule, inp::Dict{Element,<:AbstractFloat}) = inp

isunmeasured(::NullUnmeasuredRule, elm::Element) = false

"""
The ElementByStoichiometry computes the oxygen mass fraction using the standard stoichiometric rules.
"""
struct ElementByStoichiometry <: UnmeasuredElementRule 
    element::Element
    valences
    ElementByStoichiometry(elm::Element, val) = new(elm, val)
end

"""
    OByStoichiometry(valences = NeXLCore.valences)

Returns an ElementByStoichiometry for computing O using the valences provided as a tuple or Array indexed by z.
"""
OByStoichiometry(valences = NeXLCore.defaultValences) = ElementByStoichiometry(n"O", valences)

"""
    NeXLUncertainties.compute(stoic::ElementByStoichiometry, inp::Dict{Element,<:AbstractFloat})

Computes an element using stoichiometric rules.
"""
function NeXLUncertainties.compute(stoic::ElementByStoichiometry, inp::Dict{Element,<:AbstractFloat})
    inp[stoic.element] = elmbystoichiometry(stoic.element, inp, valences = stoic.valences)
    inp
end

function isunmeasured(stoic::ElementByStoichiometry, elm::Element)::Bool
    isequal(elm, stoic.element)
end

function Base.show(io::IO, stoic::ElementByStoichiometry)
    print(io,"ElementByStoichiometry[$(stoic.element)]")
end

struct ElementByFiat <: UnmeasuredElementRule
    element::Element
    massfraction::Float64
end

"""
    NeXLUncertainties.compute(stoic::ElementByFiat, inp::Dict{Element,<:AbstractFloat})

Simply assigns a quantity for an element.
"""
function NeXLUncertainties.compute(fiat::ElementByFiat, inp::Dict{Element,<:AbstractFloat})
    inp[fiat.element] = fiat.massfraction
    inp
end

function isunmeasured(fiat::ElementByFiat, elm::Element)::Bool
    isequal(elm, fiat.element)
end

function Base.show(io::IO, fiat::ElementByFiat)
    print(io,"ElementByFiat[$(fiat.element), $(fiat.massfraction)]")
end
"""
    MultiUnmeasuredElementRule

Apply multiple UnmeasuredElementRule items sequentially.
"""
struct MultiUnmeasuredElementRule <: UnmeasuredElementRule 
    rules::Vector{<:UnmeasuredElementRule}
end

function NeXLUncertainties.compute(muer::MultiUnmeasuredElementRule, inp::Dict{Element,<:AbstractFloat})
    res = inp
    for rule in muer.rules
        res = compute(rule, res)
    end
    return res
end

function isunmeasured(muer::MultiUnmeasuredElementRule, elm::Element)::Bool
    any(rule->isunmeasured(rule, elm), muer.rules)
end

function Base.show(io::IO, muer::MultiUnmeasuredElementRule)
    print(io,"MultiUnmeasuredElementRule[$(join(muer.rules, ", "))]")
end

