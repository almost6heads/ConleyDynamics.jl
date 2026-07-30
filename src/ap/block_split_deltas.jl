export is_closed_in_block, block_split_deltas, split_spectrum, connected_split_deltas

"""
    is_closed_in_block(lc::AbstractComplex, A::Vector{Int}, blockset::Set{Int})

Whether `A` is closed in the block `blockset` (Theorem 2.1: `A` is closed
in block `C` iff `cl(A) cap C = A`).
"""
function is_closed_in_block(lc::AbstractComplex, A::Vector{Int}, blockset::Set{Int})
    clA = Set(lefschetz_closure(lc, A))
    return issubset(intersect(clA, blockset), Set(A))
end

"""
    block_split_deltas(lc::AbstractComplex, block::Cells)

Enumerate every atomic split `C = A ⊔ B` of `block` with `A` closed in `C`
(Theorem 2.1, Proposition 2.2), and return the achieved
`delta = beta(A) + beta(B) - beta(C)` together with `A` and `B` for each
(the split spectrum `D({C})` of Corollary 2.3). Brute force over all
subsets, so only intended for small blocks (`length(block) <= 20`).
"""
function block_split_deltas(lc::AbstractComplex, block::Cells)
    blockI = block isa Vector{String} ? convert_cells(lc, block) : collect(block)
    n = length(blockI)
    @assert n <= 20 "block too large for brute-force split enumeration"
    blockset = Set(blockI)
    betaC = conley_index(lc, blockI)

    results = NamedTuple[]
    for mask in 1:(2^n - 2)
        bits = digits(mask, base=2, pad=n)
        A = [blockI[i] for i in 1:n if bits[i] == 1]
        if !is_closed_in_block(lc, A, blockset)
            continue
        end
        B = setdiff(blockI, A)
        betaA = conley_index(lc, A)
        betaB = conley_index(lc, B)
        delta = betaA .+ betaB .- betaC
        push!(results, (A=sort(A), B=sort(B), delta=delta))
    end
    return results
end

"""
    split_spectrum(lc::AbstractComplex, block::Cells)

The set `D({C})` of distinct delta vectors achievable by a single atomic
split of `block`.
"""
function split_spectrum(lc::AbstractComplex, block::Cells)
    return unique([r.delta for r in block_split_deltas(lc, block)])
end

"""
    connected_split_deltas(lc::AbstractComplex, block::Cells)

Like `block_split_deltas`, but restricted to splits `C = A ⊔ B` with both
`A` and `B` topologically connected (`is_connected_block`) -- the deltas
actually reachable by a single atomic refinement step while staying inside
`AP^c(X)`, as opposed to the unrestricted `AP(X)` of Theorem 5.1's
existence statement.
"""
function connected_split_deltas(lc::AbstractComplex, block::Cells)
    return filter(r -> is_connected_block(lc, r.A) && is_connected_block(lc, r.B),
                 block_split_deltas(lc, block))
end
