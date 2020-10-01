
import LinearAlgebra

"""
Define a type to hold the stability derivatives
"""
struct SDeriv
    alpha
    beta
    p
    q
    r
end
# struct SDeriv{TF}
#     alpha::TF
#     beta::TF
#     p::TF
#     q::TF
#     r::TF
# end # RMA



# overload functions using the SDeriv type
SDeriv() = SDeriv(zeros(1), zeros(1), zeros(1), zeros(1), zeros(1)) 
SDeriv(N) = SDeriv(zeros(N), zeros(N), zeros(N), zeros(N), zeros(N)) 
SDeriv(M, N) = SDeriv(zeros(M, N), zeros(M, N), zeros(M, N), zeros(M, N), zeros(M, N)) 
Base.:+(x::SDeriv, y::SDeriv) = SDeriv(x.alpha + y.alpha, x.beta + y.beta, x.p + y.p, x.q + y.q, x.r + y.r)
Base.:*(x::Any, y::SDeriv) = SDeriv(x * y.alpha, x * y.beta, x * y.p, x * y.q, x * y.r)
Base.:*(x::SDeriv, y::Any) = SDeriv(x.alpha * y, x.beta * y, x.p * y, x.q * y, x.r * y)
Base.:-(x::SDeriv) = SDeriv(-x.alpha, -x.beta, -x.p, -x.q, -x.r)
Base.:/(x::SDeriv, y::Any) = SDeriv(x.alpha / y, x.beta / y, x.p / y, x.q / y, x.r / y)
# Base.broadcast(::typeof(Base.:/), x::SDeriv, y::Any) = SDeriv(x.alpha ./ y, x.beta ./ y, x.p ./ y, x.q ./ y, x.r ./ y)
LinearAlgebra.dot(x::SDeriv, y::Array{T, 1}) where T<:Any = SDeriv(dot(x.alpha, y), dot(x.beta, y), dot(x.p, y), dot(x.q, y), dot(x.r, y))
LinearAlgebra.cross(x::SDeriv, y::Array{T, 1}) where T<:Any = SDeriv(cross(x.alpha, y), cross(x.beta, y), cross(x.p, y), cross(x.q, y), cross(x.r, y))
LinearAlgebra.cross(x::Array{T, 1}, y::SDeriv) where T<:Any = SDeriv(cross(x, y.alpha), cross(x, y.beta), cross(x, y.p), cross(x, y.q), cross(x, y.r))
Base.sum(x::SDeriv, n) = SDeriv(sum(x.alpha, n), sum(x.beta, n), sum(x.p, n), sum(x.q, n), sum(x.r, n))

function Base.getindex(x::SDeriv, i)
    return SDeriv(x.alpha[i], x.beta[i], x.p[i], x.q[i], x.r[i])
end
function Base.getindex(x::SDeriv, i, j)
    return SDeriv(x.alpha[i, j], x.beta[i, j], x.p[i, j], x.q[i, j], x.r[i, j])
end
function Base.setindex!(x::SDeriv, y::SDeriv, i)
    x.alpha[i] = y.alpha
    x.beta[i] = y.beta
    x.p[i] = y.p
    x.q[i] = y.q
    x.r[i] = y.r
end
function Base.setindex!(x::SDeriv, y::Any, i)
    x.alpha[i] = y
    x.beta[i] = y
    x.p[i] = y
    x.q[i] = y
    x.r[i] = y
end
function Base.setindex!(x::SDeriv, y::SDeriv, i, j)
    x.alpha[i, j] = y.alpha
    x.beta[i, j] = y.beta
    x.p[i, j] = y.p
    x.q[i, j] = y.q
    x.r[i, j] = y.r
end
# Base.endof(x) = length(x.alpha)