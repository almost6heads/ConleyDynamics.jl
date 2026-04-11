export lefschetz_gfp_conversion

"""
    lcgfp = lefschetz_gfp_conversion(lc::AbstractComplex, p::Int)

Convert a Lefschetz complex to the same complex over a finite field.

It is expected that the boundary matrix of the given Lefschetz
complex `lc` is defined over the rationals, and that the target
characteristic `p` is a prime.
"""
function lefschetz_gfp_conversion(lc::AbstractComplex, p::Int)
    #
    # Convert a Lefschetz complex to a finite field complex
    #

    # Extract the boundary matrix data

    nr, nc, tchar, tzero, tone, r, c, v = lists_from_sparse(lc.boundary)

    # Make sure all the requirements are met

    if !(tone isa Rational{Int})
        error("The Lefschetz complex has to be over the rationals!")
    end

    if p <= 0
        error("The target characteristic has to be a positive prime!")
    end

    # Make sure that the boundary matrix entries are in fact integers

    if maximum(denominator.(v)) > 1
        error("The boundary matrix contains fractional entries!")
    end

    # Convert the entries to integers and create the new boundary matrix

    vint = convert(Vector{Int},v)
    boundarygfp = sparse_from_lists(nr,nc,p,Int(0),Int(1),r,c,vint)

    # Create the new Lefschetz complex and return it

    lcgfp = LefschetzComplex(lc.labels, lc.dimensions, boundarygfp; validate=false)
    return lcgfp
end

