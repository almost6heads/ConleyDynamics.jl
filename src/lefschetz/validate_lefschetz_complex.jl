export validate_lefschetz_complex

"""
    validate_lefschetz_complex(lc::AbstractComplex)

Check whether the boundary matrix of `lc` satisfies `∂² = 0`.

Returns `true` if the complex is valid. Throws an error if the boundary
matrix does not square to zero, indicating an incorrectly assembled complex.

This function is the explicit counterpart to the `validate=true` keyword
argument of the [`LefschetzComplex`](@ref) constructor. It is useful for
checking a complex that was not validated at construction time.

# Examples
```jldoctest
julia> lc = create_simplicial_complex(["A","B","C"], [["A","B"],["B","C"]]);

julia> validate_lefschetz_complex(lc)
true
```
"""
function validate_lefschetz_complex(lc::AbstractComplex)
    #
    # Check whether the boundary matrix squares to zero
    #
    if !sparse_is_zero(sparse_multiply(lc.boundary, lc.boundary))
        error("Boundary matrix does not square to zero: invalid Lefschetz complex.")
    end
    return true
end
