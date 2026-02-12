export example_torsion_chaos

"""
    sc, vfG, vfC = example_torsion_chaos(n::Int, p::Int)

Create a triangulation of a space with 1-dimensional torsion, as well
as two Forman vectr fields on this complex.

The function returns a simplicial complex `sc` which has the following
integer homology groups:

- In dimension 0 it is the group of integers.
- In dimension 1 it is the integers modulo `n`.
- In dimension 2 it is the trivial group.

In other words, the simplicial complex has nontrivial torsion
in dimension 1. It is a triangulation of an `n`-gon, in which
all boundary edges are oriented counterclockwise, and all of
these edges are identified. The parameter `p` specifies the
characteristic of the underlying field.

In addition, two Forman vector fields `vfG` and `vfC` are returned.
The first one is a gradient vector field whose connection matrix
has a large connection matrix entry. In fact, if `p` is any prime
larger than `n` then there will be an entry `n` in the matrix.
The scond Forman vector field contains a chaotic Morse set. This
Morse set will have trivial Morse index for most `p`. On the other
hand, for prime `p = n` the set has the Morse index of an unstable 
periodic orbit.

# Examples
```jldoctest
julia> sc, vfG, vfC = example_torsion_chaos(3,7);

julia> homology(sc)
3-element Vector{Int64}:
 1
 0
 0

julia> cmG = connection_matrix(sc, vfG);

julia> sparse_show(cmG.matrix)
 0 0 0
 0 0 3
 0 0 0

julia> print(cmG.labels)
["0w", "0w0x", "0w0x1s"]

julia> cmC = connection_matrix(sc, vfC);

julia> sparse_show(cmC.matrix)
 0

julia> print(cmC.labels)
["0w"]

julia> msC, psC = morse_sets(sc, vfC, poset=true);

julia> [conley_index(sc, mset) for mset in msC]
2-element Vector{Vector{Int64}}:
 [1, 0, 0]
 [0, 0, 0]

julia> psC
2×2 Matrix{Bool}:
 0  1
 0  0

julia> length.(msC)
2-element Vector{Int64}:
  1
 26
```
"""
function example_torsion_chaos(n::Int, p::Int)
    # 
    # Create a triangulation of a space with 1-d torsion,
    # as well as two Forman vector fields on it.
    #
    
    # Create the simplicial complex

    sc = simplicial_torsion_space(n,p)

    # Compute the number of digits for the labels

    labwidth = Int(ceil(log(0.5+n)/log(10)))
    labformat = "%0" * string(labwidth) * "d%s"

    # Create the gradient vector field

    vfG = Vector{Vector{String}}()

    # First the arrows in the center of the triangulation

    vlabelc = Printf.format(Printf.Format(labformat), 0, "c")
    for k = 0:n-1
        vlabels  = Printf.format(Printf.Format(labformat), k+1, "s")
        vlabelt  = Printf.format(Printf.Format(labformat), k+1, "t")
        kn = mod(k+1,n)
        vlabelsn = Printf.format(Printf.Format(labformat), kn+1, "s")

        if k==0
            push!(vfG, [vlabelc, join(sort([vlabelc, vlabels]))])
        else
            push!(vfG, [vlabels, join(sort([vlabelc, vlabels]))])
        end
        push!(vfG, [vlabelt, join(sort([vlabelc, vlabelt]))])
        push!(vfG, [join(sort([vlabels,  vlabelt])), join(sort([vlabels,  vlabelt, vlabelc]))])
        push!(vfG, [join(sort([vlabelsn, vlabelt])), join(sort([vlabelsn, vlabelt, vlabelc]))])
    end

    # Now the arrows around the outer layer

    vlabelb1 = Printf.format(Printf.Format(labformat), 0, "w")
    vlabelb2 = Printf.format(Printf.Format(labformat), 0, "x")
    vlabelb3 = Printf.format(Printf.Format(labformat), 0, "y")
    vlabelb4 = Printf.format(Printf.Format(labformat), 0, "z")
    for k = 0:n-1
        vlabels  = Printf.format(Printf.Format(labformat), k+1, "s")
        vlabelt  = Printf.format(Printf.Format(labformat), k+1, "t")
        kn = mod(k+1,n)
        vlabelsn = Printf.format(Printf.Format(labformat), kn+1, "s")

        if k==0
            push!(vfG, [vlabels, join(sort([vlabels, vlabelb1]))])
        else
            push!(vfG, [join(sort([vlabels, vlabelb1])), join(sort([vlabels, vlabelb1, vlabelb2]))])
        end
        push!(vfG, [join(sort([vlabels,  vlabelb2])), join(sort([vlabels,  vlabelb2, vlabelb3]))])
        push!(vfG, [join(sort([vlabels,  vlabelb3])), join(sort([vlabels,  vlabelt,  vlabelb3]))])
        push!(vfG, [join(sort([vlabelt,  vlabelb3])), join(sort([vlabelt,  vlabelb3, vlabelb4]))])
        push!(vfG, [join(sort([vlabelt,  vlabelb4])), join(sort([vlabelt,  vlabelsn, vlabelb4]))])
        push!(vfG, [join(sort([vlabelsn, vlabelb4])), join(sort([vlabelsn, vlabelb1, vlabelb4]))])
    end

    # Next the arrows around the repeated outer edges

    push!(vfG, [vlabelb2, vlabelb2 * vlabelb3])
    push!(vfG, [vlabelb3, vlabelb3 * vlabelb4])
    push!(vfG, [vlabelb4, vlabelb1 * vlabelb4])

    # Finally, create the chaotic Forman vector field

    vfC = deepcopy(vfG)
    vlabels = Printf.format(Printf.Format(labformat), 1, "s")
    push!(vfC, [vlabelb1 * vlabelb2, vlabelb1 * vlabelb2 * vlabels])

    # Return the results

    return sc, vfG, vfC
end

