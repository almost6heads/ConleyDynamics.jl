export lefschetz_clomo_pair

"""
    lefschetz_clomopair(lc::LefschetzComplex, subcomp::Vector{Int})

Determine the closure-mouth-pair associated with a Lefschetz complex subset.

The function returns the pair `(closure,mouth)`.
"""
function lefschetz_clomo_pair(lc::LefschetzComplex, subcomp::Vector{Int})
    #
    # Determine the closure-mouth-pair for a Lefschetz complex subset
    #

    subclosure = lefschetz_closure(lc, subcomp)

    return subclosure, setdiff(subclosure, subcomp)
end

"""
    lefschetz_clomopair(lc::LefschetzComplex, subcomp::Vector{String})

Determine the closure-mouth-pair associated with a Lefschetz complex subset.

The function returns the pair `(closure,mouth)`.
"""
function lefschetz_clomo_pair(lc::LefschetzComplex, subcomp::Vector{String})
    #
    # Determine the closure-mouth-pair for a Lefschetz complex subset
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end
    subclosureI, submouthI = lefschetz_clomo_pair(lc, subcompI)
    subclosure = lc.labels[subclosureI]
    submouth   = lc.labels[submouthI]

    return subclosure, submouth
end

