"""
    const MakieType = UnionMakieType

Deprecated parent type for dispatch on Makie plotters.
"""
const MakieType = UnionMakieType
Base.@deprecate_binding MakieType UnionMakieType

export MakieType
