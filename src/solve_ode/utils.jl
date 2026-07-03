"""
    _muladd!(A, cst, B, C)

Compute in-place `A = cst * B + C` by iterating directly over the nonzero values.
Assumes `A`, `B`, and `C` share the same sparsity pattern.
"""
function _muladd!(
        A::Symmetric{T, SparseMatrixCSC{T, I}},
        cst::T,
        B::Symmetric{T, SparseMatrixCSC{T, I}},
        C::Symmetric{T, SparseMatrixCSC{T, I}}) where {T, I}
    nzval_A = A.data.nzval
    nzval_B = B.data.nzval
    nzval_C = C.data.nzval

    @. nzval_A = muladd(cst, nzval_B, nzval_C)

    return nothing
end

function update_state!(
        state::State{T},
        v̂ⁿ::AbstractArray{T},
        ĉⁿ::AbstractArray{T},
        r̂ⁿ::AbstractArray{T},
        τ::T
) where {T}
    state.n += 1
    state.t += τ

    @. state.v = muladd(2, v̂ⁿ, -state.v)
    @. state.d = muladd(τ, v̂ⁿ, state.d)

    @. state.c = muladd(2, ĉⁿ, -state.c)

    @. state.r = muladd(2, r̂ⁿ, -state.r)
    @. state.z = muladd(τ, r̂ⁿ, state.z)

    return nothing
end