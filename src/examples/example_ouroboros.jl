export example_ouroboros

"""
    ec, vf = example_ouroboros()

Create the Ouroboros example of Justin Thorpe.

The example consists of a triangulation `ec` of the rectangle
`[0,4] x [0,5]`, which is created by `create_simplicial_rectangle(4, 5,
euclidean=true)`, and a planar vector field `vf`. The vector field points
toward the point `(2.1, 2.5)`, except near a few edge midpoints, where it
is modified so that the flow-induced multivector field contains a
multivector whose closure winds around the center and touches itself at a
single vertex. The closure of this multivector is an isolating block whose
entry set is not closed, i.e., it is not a strong isolating block. This
example appears in the dissertation *Periodic Orbits in the Setting of
Combinatorial Dynamics* by Justin Thorpe (George Mason University, 2026),
and it was contributed in pull request #10.

The flow enters the rectangle through every boundary edge. The
accumulating multivector field contains the Ouroboros multivector, while
the finest strongly accumulating multivector field merges it with one of
its neighbors.

# Example

```julia
ec, vf = example_ouroboros()
fwd = planar_transitions(ec, vf)
bwd = planar_transitions(ec, x -> -vf(x))
lc  = euclidean_to_lefschetz(ec)
mvf = create_mvf_hull(lc, fwd)
duals = mvf_dual(lc, mvf, bwd)
```
"""
function example_ouroboros()
    #
    # Create the Ouroboros example
    #
    ec = create_simplicial_rectangle(4, 5, euclidean=true)
    return ec, ouroboros_vf
end

function ouroboros_vf(x::Vector{Float64})
    #
    # The vector field of the Ouroboros example (by Justin Thorpe)
    #
    if norm(x - [2.0, 3.5]) < 0.1
        return [0.0, 1.0]
    elseif norm(x - [2.25, 3.25]) < 0.1
        return [1.0, 1.0]
    elseif norm(x - [2.5, 3.0]) < 0.1
        return [1.0, 0.0]
    elseif norm(x - [2.75, 2.75]) < 0.1
        return [1.0, 1.0]
    elseif norm(x - [2.75, 2.25]) < 0.1
        return [1.0,-1.0]
    elseif norm(x - [2.5, 2.0]) < 0.1
        return [1.0, 0.0]
    elseif norm(x - [2.25, 1.75]) < 0.1
        return [1.0,-1.0]
    elseif norm(x - [2.0, 1.5]) < 0.1
        return [0.0,-1.0]
    elseif norm(x - [1.75, 1.75]) < 0.1
        return [1.0, 1.0]
    elseif norm(x - [1.5, 2.0]) < 0.1
        return [1.0, 0.0]
    elseif norm(x - [1.25, 2.25]) < 0.1
        return [1.0, 1.0]
    elseif norm(x - [1.25, 2.75]) < 0.1
        return [1.0,-1.0]
    elseif norm(x - [0.25, 2.25]) < 0.1
        return [1.0, 1.0]
    elseif norm(x - [0.75, 2.25]) < 0.1
        return [1.0,-1.0]
    else
        return [2.1, 2.5] - x
    end
end
