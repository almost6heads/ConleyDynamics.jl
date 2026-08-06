export morse_vector, block_rho

"""
    _rank_decomposition(d::Vector{Int})

Solve `d_k = rho_k + rho_{k+1}` for `rho`, given `rho_0 = 0` and `d`
indexed `1:n+1` for dimensions `0:n`. The same recursion also extracts
the `c_k` in the canonical decomposition of a
Morse-vector jump `delta = sum c_k(e_k+e_{k-1})`.

Internal helper, not exported.
"""
function _rank_decomposition(d::Vector{Int})
    n   = length(d) - 1
    rho = zeros(Int, n + 1)
    for k in 1:n
        rho[k+1] = d[k] - rho[k]
    end
    return rho
end

"""
    morse_vector(lc::AbstractComplex, mvf::CellSubsets)

Morse vector `M(W) = sum_{V in W} beta(V)` of a multivector field/partition
`mvf`. Every non-critical block contributes 0, so this is just the sum of
the Conley indices of the critical multivectors.
"""
function morse_vector(lc::AbstractComplex, mvf::CellSubsets)
    crit = mvf_critical(lc, mvf)
    M = zeros(Int, lc.dim + 1)
    for mv in crit
        M .+= conley_index(lc, mv)
    end
    return M
end

"""
    block_rho(lc::AbstractComplex, block::Cells)

The rank vector `(rho_0(V),...,rho_n(V))` of a single locally closed
multivector `block`, recovered from `beta(V)` and `f(V)`:
`f(V) - beta(V) = sum_k rho_k(V) (e_k + e_{k-1})`.

Uses `beta(V) = conley_index(lc, block)`, which coincides with the
absolute homology of the restricted boundary operator for any locally
closed `V`.
"""
function block_rho(lc::AbstractComplex, block::Cells)
    blockI = block isa Vector{String} ? convert_cells(lc, block) : collect(block)
    n = lc.dim
    f = zeros(Int, n + 1)
    for c in blockI
        f[lc.dimensions[c]+1] += 1
    end
    beta = conley_index(lc, blockI)
    d    = f .- beta
    rho  = _rank_decomposition(d)
    @assert rho[end] == d[end] "Inconsistent rank recursion for block $blockI"
    return rho
end
