# --- reference quantities --- #

"""
    Reference(S, c, b, r)

Reference quantities.

**Arguments**
 - `S`: reference area
 - `c`: reference chord
 - `b`: reference span
 - `r`: reference location for all rotations/moments
 - `V`: reference velocity (magnitude)
 - `R`: basis vectors of the vehicle frame expressed in global coordinates
"""
struct Reference{TF}
    S::TF
    c::TF
    b::TF
    r::SVector{3, TF}
    V::TF
end

function Reference(S, c, b, r, V)
    TF = promote_type(typeof(S), typeof(c), typeof(b), eltype(r), typeof(V))
    return Reference{TF}(S, c, b, r, V)
end

Base.eltype(::Type{Reference{TF}}) where TF = TF
Base.eltype(::Reference{TF}) where TF = TF

Reference{TF}(r::Reference) where TF = Reference{TF}(r.S, r.c, r.b, r.r, r.V)
Base.convert(::Type{Reference{TF}}, r::Reference) where {TF} = Reference{TF}(r)

# --- reference frames --- #

"""
    AbstractFrame

Supertype for the different possible reference frames used by this package.
"""
abstract type AbstractFrame end

"""
   Body <: AbstractFrame

Reference frame aligned with the global X-Y-Z axes
"""
struct Body <: AbstractFrame end

"""
    Stability <: AbstractFrame

Reference frame rotated from the body frame about the y-axis to be aligned with
the freestream `alpha`.
"""
struct Stability <: AbstractFrame end

"""
    Wind <: AbstractFrame

Reference frame rotated to be aligned with the freestream `alpha` and `beta`
"""
struct Wind <: AbstractFrame end
