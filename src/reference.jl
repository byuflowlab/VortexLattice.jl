"""
    Reference(S, c, b, rcg)

Reference quantities.

**Arguments**
- `S`: reference area
- `c`: reference chord
- `b`: reference span
- `rcg`: location of center of gravity
"""
struct Reference{TF}
    S::TF
    c::TF
    b::TF
    rcg::SVector{3, TF}
end

function Reference(S, c, b, rcg)
    TF = promote_type(typeof(S), typeof(c), typeof(b), eltype(rcg))
    return Reference{TF}(S, c, b, rcg)
end

Base.eltype(::Type{Reference{TF}}) where TF = TF
Base.eltype(::Reference{TF}) where TF = TF
