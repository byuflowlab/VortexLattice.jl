"""
    Reference(S, c, b, r)

Reference quantities.

**Arguments**
- `S`: reference area
- `c`: reference chord
- `b`: reference span
- `r`: reference location for all rotations/moments
"""
struct Reference{TF}
    S::TF
    c::TF
    b::TF
    r::SVector{3, TF}
end

function Reference(S, c, b, r)
    TF = promote_type(typeof(S), typeof(c), typeof(b), eltype(r))
    return Reference{TF}(S, c, b, r)
end

Base.eltype(::Type{Reference{TF}}) where TF = TF
Base.eltype(::Reference{TF}) where TF = TF
