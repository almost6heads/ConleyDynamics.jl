
# Tip: Release Notes

Did you know you can add release notes too? Just add markdown formatted text
underneath the comment after the text "Release notes:" and it will be added to
the registry PR, and if TagBot is installed it will also be added to the release
that TagBot creates. i.e., you could add a "## Breaking Changes" header with a
description, etc.

*******************************************************************************

@JuliaRegistrator register

Release notes:

## v0.3.1 (June 7, 2025)

- Added `lefschetz_reduction_maps`. It computes the reduction of a Lefschetz
  complex based on a sequence of elementary reduction pairs, but also provides
  the involved chain equivalences and chain homotopies.
- Added the new sparse matrix function `sparse_zero`, and changed the name
  of `sparse_nonzero_count` to `sparse_nz_count`.
- Modified `sparse_set_entry` in the finite field case.
- Added the helper function `scalar_inverse`.

## v0.3.0 (June 6, 2025)

This release does not contain any breaking changes. But the 
following new Lefschetz complex functionality and examples
have been added since release 0.2.0:

- Added `lefschetz_interior` to determine the interior of a Lefschetz
  complex subset, if the Lefschetz complex is interpreted as a finite
  topological space.
- Added `lefschetz_topboundary` to determine the topological boundary
  in the above setting.
- Added `example_dunce_chaos`. It constructs a Forman vector field
  on a minimal representation of the dunce hat which has a chaotic
  isolated invariant set with trivial Conley index.
- Added `example_torsion_chaos`. It constructs a Forman vector field
  on a simplicial complex with torsion which has a chaotic isolated
  invariant set with trivial Conley index, for certain values of the
  field characteristic. In addition, an associated gradient vector
  field has large entries in the connection matrix.

In addition, and in anticipation of future extensions, more sparse
matrix functions have been added:

- The functions `sparse_add`, `sparse_subtract`, and `sparse_scale`
  perform sparse matrix addition, subtraction, and scalar multiplication.
  All of these functions can also be involved using the operators
  `+`, `-`, and `*`, respectively.
- A new function `sparse_nonzero_count` determines the number of
  nonzero entries of a sparse matrix.

## v0.2.4 (June 2, 2025)

- Added documentation for the functions `example_dunce_chaos`
  and `example_torsion_chaos`.

## v0.2.3 (May 12, 2025)

- Added `example_torsion_chaos`. It constructs a Forman vector field
  on a simplicial complex with torsion which has a chaotic isolated
  invariant set with trivial Conley index, for certain values of the
  field characteristic. In addition, an associated gradient vector
  field has large entries in the connection matrix.
- The documentation for this function still has to be written.

## v0.2.2 (May 11, 2025)

- Added `example_dunce_chaos`. It constructs a Forman vector field
  on a minimal representation of the dunce hat which has a chaotic
  isolated invariant set with trivial Conley index.
- The documentation for this function still has to be written.

## v0.2.1 (March 27, 2025)

- Added `lefschetz_interior` to determine the interior of a Lefschetz
  complex subset, if the Lefschetz complex is interpreted as a finite
  topological space.
- Added `lefschetz_topboundary` to determine the topological boundary
  in the above setting.

## v0.2.0 (March 7, 2025)

This release does not contain any breaking changes. But the 
following new functionality has been added since release 0.1.0:

- A function to easily create small Lefschetz complexes over `GF(2)`.
- Reduction of Lefschetz complexes via elementary reductions.
- Basic functionality to work with filters and shallow pairs.
- New functions to create triangulations for a number of
  sample topological spaces. These are mostly for demonstration
  purposes during teaching.

## v0.1.7 (March 4, 2025)

- Added `lefschetz_information` for basic information about a
  Lefschetz complex.
- Updated the documentation.

## v0.1.6 (March 1, 2025)

- Added triangulations for the torus, the Klein bottle and the
  projective plane.
- Updated the documentation.

## v0.1.5 (February 25, 2025)

- Added `create_random_filter` to create a random injective filter.
- Added `filter_shallow_pairs` to find all shallow pairs of a filter.
- Added `filter_induced_mvf` to compute the multivector field induced
  by a filter.

## v0.1.4 (February 21, 2025)

- Added `lefschetz_reduction` to construct a smaller Lefschetz complex
  via elementary reductions.
- Added `sparse_get_nz_row` to extract all column indices corresponding
  to nonzero entries in a row.
- Added `sparse_is_zero` to check whether a sparse matrix is zero.
- Added consistency check to `create_lefschetz_gf2`.

## v0.1.3 (February 17, 2025)

- Added `create_lefschetz_gf2` function for quick Lefschetz complex creation
  over a field of characteristic 2.
- Simplified some examples in the documentation based on the new function.

## v0.1.2 (February 12, 2025)

- Fixed a rare exception in the multivector field creation

## v0.1.1 (December 7, 2024)

- Removed unnecessary test from Lefschetz complex creation for speedup reasons.
- Minor changes to the documentation.

## v0.1.0 (November 29, 2024)

This is the initial official version hosted on the general Julia registry.

