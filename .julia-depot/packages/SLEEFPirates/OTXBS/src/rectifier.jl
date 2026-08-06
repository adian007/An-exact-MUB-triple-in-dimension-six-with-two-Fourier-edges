
struct PReLu{T}
  a::T
end
@inline (p::PReLu)(x) = ifelse(x <= 0, x * p.a, x)
# Convert integer inputs to the matching float type up front so the kernel
# runs in a consistent precision. Without this the integer path would mix
# `Int` with a `Float64` produced by `Base.exp(::Integer)` and then narrow
# at the end, giving a 1-ULP delta vs the pure-float path. See `Φ` below
# for the same pattern.
@inline (p::PReLu)(x::IntegerType) = p(float(x))

@inline Φ(x::T) where {T<:FloatType} =
  @fastmath T(0.5) * (one(T) + VectorizationBase.verf(x * T(0.7071067811865476)))
@inline Φ(x::IntegerType) = Φ(float(x))
@inline gelu(x) = Base.FastMath.mul_fast(x, Φ(x))
@inline gelu(x::IntegerType) = gelu(float(x))

@inline softplus(x) = log1p(Base.exp(x))
@inline softplus(x::IntegerType) = softplus(float(x))
@inline silu(x) = x * sigmoid_fast(x)
@inline silu(x::IntegerType) = silu(float(x))

struct Elu{T}
  a::T
end
@inline (e::Elu)(x) = ifelse(x <= 0, Base.FastMath.mul_fast(e.a, expm1_fast(x)), x)
@inline (e::Elu)(x::IntegerType) = e(float(x))
