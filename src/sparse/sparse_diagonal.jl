export sparse_diagonal

"""
    sparse_diagonal(nr::Int, nc::Int, d::Vector{<:Real};
                    p::Int=0, offset::Int=0)

Create sparse diagonal matrix.

The function creates a sparse matrix with `nr` rows and `nc`
columns, over a field with characteristic `p`. Without additional
arguments, the matrix has the entries in the vector `d` along the
main diagonal. If the optional argument `offset` is given and
positive, then the diagonal is placed `offset` positions above
the main diagonal, if it is negative, then it indicates how many
positions below the diagonal it is placed. The function places as
many entries from `d` as possible, i.e., if the length of `d` is
too short, the remaining diagonal entries are zero, if it is too
long, then the later entries are ignored. The function tries to
convert the entries in `d` to the correct format based on `p`.
Errors are raised if this is not possible.
"""
function sparse_diagonal(nr::Int, nc::Int, d::Vector{<:Real};
                         p::Int=0, offset::Int=0)
    #
    # Create a sparse diagonal matrix
    #

    # Initialize some variables

    tchar = p
    if tchar == 0
        tone  = Rational(1)
        tzero = Rational(0)
    else
        tone  = Int(1)
        tzero = Int(0)
    end

    # Create the index and entry lists

    r = Vector{Int}()
    c = Vector{Int}()
    v = Vector{typeof(tzero)}()

    # Initialize the indices

    if offset >= 0
        ri = 1
        ci = 1 + offset
    else
        ri = 1 - offset
        ci = 1
    end

    # Add the diagonal entries

    for dentry in d
        if (1 <= ri <= nr) && (1 <= ci <= nc)
            if tchar == 0
                push!(r, ri)
                push!(c, ci)
                push!(v, Rational(dentry))
            else
                push!(r, ri)
                push!(c, ci)
                push!(v, Int(dentry))
            end
        end
        ri += 1
        ci += 1
    end

    # Return the result

    return sparse_from_lists(nr, nc, tchar, tzero, tone, r, c, v)
end

