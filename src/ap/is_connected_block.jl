export is_connected_block, is_connected_partition

"""
    is_connected_block(lc::AbstractComplex, block::Cells)

Whether `block` is connected as a subspace of `lc` -- the "c" in `AP^c(X)`.

This is a genuinely different question from `block_rho`'s `rho_0`/`beta_0`:
those come from the *restricted boundary operator* on the block (Fact 1.1)
and measure homological cancellation, which is 0 for any regular
(non-critical) multivector regardless of whether it is topologically
connected -- e.g. a Forman pair `{v,e}` with `v` a face of `e` is connected
as a subspace but has `beta(V) = (0,0,...)` because it is non-critical.
Connectivity here instead means: is the comparability graph of `block`
(edges wherever one cell is a direct face of another, both cells inside
`block`) a single connected component. Reachability via direct face
relations already captures full order-comparability, since any longer
chain `x < y < z` decomposes into direct-face steps.
"""
function is_connected_block(lc::AbstractComplex, block::Cells)
    blockI = block isa Vector{String} ? convert_cells(lc, block) : collect(block)
    n = length(blockI)
    n <= 1 && return true

    blockset = Set(blockI)
    adj = Dict(c => Int[] for c in blockI)
    for c in blockI
        for f in lefschetz_boundary(lc, c)
            if f in blockset
                push!(adj[c], f)
                push!(adj[f], c)
            end
        end
    end

    visited = Set{Int}([blockI[1]])
    stack   = [blockI[1]]
    while !isempty(stack)
        x = pop!(stack)
        for y in adj[x]
            if !(y in visited)
                push!(visited, y)
                push!(stack, y)
            end
        end
    end

    return length(visited) == n
end

"""
    is_connected_partition(lc::AbstractComplex, mvf::CellSubsets)

Whether every explicit multivector of `mvf` is connected
(`is_connected_block`); implicit singletons are always connected. This is
exactly the membership test for `AP^c(X)`, the sub-poset that
`construct_ap_space(lc; connected=true)` enumerates directly.
"""
function is_connected_partition(lc::AbstractComplex, mvf::CellSubsets)
    return all(b -> is_connected_block(lc, b), mvf)
end
