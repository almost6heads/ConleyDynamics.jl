export example_moebius

"""
    lc1, mvf1, lc2, mvf2 = example_moebius()

Create two simplicial complexes for a cylinder and Moebius
strip, respectively, together with associated multivector
fields on them.

The multivector field is the same, and it has one critical 
cell each in dimension 1 and 2 in the interior of the strip.
The boundary consists of two periodic orbits for `lc1` and
`mvf1`, and of one periodic orbit in the Moebius case `lc2`
and `mvf2`. The latter case leads to different connection
matrices for the fields GFP(2) and GFP(7), for example.

# Examples
```jldoctest
julia> lc1, mvf1, lc2, mvf2 = example_moebius();

julia> cmp2 = connection_matrix(lc2, mvf2; p=2);

julia> cmp7 = connection_matrix(lc2, mvf2; p=7);

julia> cmp2.cm
[0   0   0   0]
[0   0   0   1]
[0   0   0   0]
[0   0   0   0]

julia> cmp7.cm
[0   0   0   0]
[0   0   0   1]
[0   0   0   2]
[0   0   0   0]
```
"""
function example_moebius()
    #
    # Example which returns different connection matrices for
    # different finite fields. The two examples are a cylinder
    # and a Moebius strip. The second example has connection 
    # matrices which differ for p=2 and p=7.
    #

    # Create the labels and simplices

    labels = ["A","B","C","D","E","F","G","H"]
    simplices1 = [["A","B","C"],["B","C","D"],["C","D","E"],["D","E","F"],
                  ["E","F","G"],["F","G","H"],["A","G","H"],["A","B","H"]]
    simplices2 = [["A","B","C"],["B","C","D"],["C","D","E"],["D","E","F"],
                  ["E","F","G"],["F","G","H"],["B","G","H"],["A","B","H"]]

    # Define the multivector fields. They are the same, except for the
    # twisting of the underlying simplicial complex in the second example.

    mvf1 = [["A","AC"],["C","CE"],["E","EG"],["G","AG"],
            ["B","BD"],["D","DF"],["F","FH"],["H","BH"],
            ["EF","DEF"],["DE","CDE"],["CD","BCD"],["BC","ABC"],
            ["FG","FGH"],["GH","AGH"],["AH","ABH"]]
    mvf2 = [["A","AC"],["C","CE"],["E","EG"],["G","BG"],
            ["B","BD"],["D","DF"],["F","FH"],["H","AH"],
            ["EF","DEF"],["DE","CDE"],["CD","BCD"],["BC","ABC"],
            ["FG","FGH"],["GH","BGH"],["BH","ABH"]]

    # Create the simplicial complexes

    lc1 = create_simplicial_complex(labels, simplices1)
    lc2 = create_simplicial_complex(labels, simplices2)

    # Return the complexes and multivector fields

    return lc1, mvf1, lc2, mvf2
end
