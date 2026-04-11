export lefschetz_field

"""
    fieldstr = lefschetz_field(lc::AbstractComplex)

Returns the Lefschetz complex coefficient field.
"""
function lefschetz_field(lc::AbstractComplex)
    #
    # Return the Lefschetz complex coefficient field
    #

    p = lc.boundary.char

    if (p == 0) && (lc.boundary.one isa Rational{Int})
        fieldstr = "QQ"
    elseif (p > 0) && (lc.boundary.one isa Int)
        fieldstr = "GF(" * string(p) * ")"
    else
        fieldstr = "unknown!"
    end

    return fieldstr
end

