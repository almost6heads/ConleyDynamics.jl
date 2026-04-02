
# Tip: Release Notes

Did you know you can add release notes too? Just add markdown formatted text
underneath the comment after the text "Release notes:" and it will be added to
the registry PR, and if TagBot is installed it will also be added to the release
that TagBot creates. i.e., you could add a "## Breaking Changes" header with a
description, etc.

*******************************************************************************

@JuliaRegistrator register

Release notes:

## v0.6.3 (April 2, 2026)

- Added the function `mvf_is_acyclic`.
- Added the function `lefschetz_neighbors` which determines
  all cells that can be added to a cell subset without
  increasing its number of connected components.

## v0.6.2 (March 31, 2026)

- Added the functions `sparse_basis_kernel` and `sparse_basis_range`
  to find bases for kernel and range of a sparse matrix.
- Added `sparse_solve` to find a particular solution of a sparse
  linear system.
- Added the test `sparse_is_rref`.

## v0.6.1 (March 30, 2026)

- Added the functions `sparse_rref` and `sparse_rref!` to
  compute the reduced row Echelon form of a sparse matrix.
- Added a method to `sparse_from_full`, which allows for
  an integer vector argument and creates a sparse matrix
  in column form.
- Added the function `sparse_size` with only one imatrix
  argument, which returns the number of rows and columns
  of the matrix as a pair.

## v0.6.0 (March 24, 2026)

This release does not contain any breaking changes. But the 
following functionality has been added since release 0.5.0:

- There are now four different algorithms implemented for the
  computation of connection matrices.
- There are now three different algorithms implemented for the
  computation of persistence.
- New functions for the treatment of Forman gradient vector
  fields such as `forman_conley_maps`, `forman_critical_cells`,
  and `forman_all_cell_types`.
- New functions for sparse matrices include concatenation
  functions, `sparse_diagonal`, and `sparse_transpose`.
  In addition, the output format based on `sparse_show` has
  been improved.
- A number of functions have been rewritten for speed, and 
  parallelized whenever possible.

## v0.5.8 (March 22, 2026)

- Added the function `sparse_diagonal`.
- Reorganized the sparse matrix section of the manual and for the API.

## v0.5.7 (March 21, 2026)

- Added the function `sparse_hvcat`.
- Added the function `sparse_hvncat`.
- Added the function `sparse_cat`.
- One can now use Julia's block matrix notation
  to form large sparse matrices using smaller blocks.

## v0.5.6 (March 19, 2026)

- Added the function `sparse_vcat`.
- Added new methods to `Base.hcat` and `Base.vcat` so that the 
  concatenation functions can use the usual short array syntax.

## v0.5.5 (March 19, 2026)

- Added the function `sparse_transpose`.
- Added the function `sparse_hcat`.

## v0.5.4 (March 11, 2026)

- Changed the output format for sparse matrices.
- Added the function `forman_conley_maps`.
- Added the helper functions `forman_critical_cells` and `forman_all_cell_types`.

## v0.5.3 (March 5, 2026)

- Some under-the-hood speed improvements.
- There are now three different algorithms implemented
  for the computation of homology and persistent homology.
- There are an additional two algorithms for the
  computation of connection matrices.

## v0.5.2 (March 4, 2026)

Changed the algorithm default.

## v0.5.1 (March 3, 2026)

Some quick additions.. :-)

## v0.5.0 (February 28, 2026)

This release does not contain any breaking changes. But the 
following functionality has been added since release 0.4.0:

- The connection matrix can now be computed with two different
  algorithms. In addition to the original `DLMS24` algorithm,
  one can now also use the faster `DHL26` algorithm. The latter
  one is now the default.
- Added the function `forman_gpaths` to find all gradient paths in
  a Forman gradient vector field.
- Added `forman_path_weight`, which computes the weight of a Forman
  gradient path.
- Added `lefschetz_filtration_mvf` to determine the multivector field
  associated with a filtration on a Lefschetz complex.
- Added `cellsubset_distance` to determine the distance from a cellsubset
  to a given point, provided coordinates for the vertices are provided.
- Added `cellsubset_planar_area` to compute the area of a cellsubset
  in a planar Lefschetz complex. This function does allow for 2-cells
  that are polygonal.

## v0.4.4 (February 26, 2026)

Bug fix release..

## v0.4.3 (February 26, 2026)

- Added the function `forman_gpaths` to find all gradient paths in
  a Forman gradient vector field.
- Added `forman_path_weight`, which computes the weight of a Forman
  gradient path.
- Added the helper functions `scalar_multiply` and `scalar_add` to allow
  for easier computations in the underlying field.

## v0.4.2 (February 19, 2026)

- Added `lefschetz_filtration_mvf` to determine the multivector field
  associated with a filtration on a Lefschetz complex.
- Added `cellsubset_distance` to determine the distance from a cellsubset
  to a given point, provided coordinates for the vertices are provided.
- Added `cellsubset_planar_area` to compute the area of a cellsubset
  in a planar Lefschetz complex. This function does allow for 2-cells
  that are polygonal.

## v0.4.1 (February 16, 2026)

Maintenance release: Unified float argument type declarations

## v0.4.0 (February 15, 2026)

This release does not contain any breaking changes. But the 
following functionality has been added since release 0.3.0:

- New functions for selecting cell subsets from a list
  of cell subsets, based on their location with respect to
  a rectangle or circle in the plane.
- Functions which allow for graded basis changes in Lefschetz 
  complexes.
- Functions which implement reduction pairs for Lefschetz
  complexes. These also provide the resulting chain maps
  and chain homotopies.
- Functions for computing Forman's combinatorial flow and
  its stabilization.
- Functions for generating and extracting chains for the
  underlying Lefschetz chain complexes.
- New sparse matrix functionality, such as comparison,
  computing the inverse, and printing with row and 
  column labels.

## v0.3.11 (February 13, 2026)

- Added `Base.show` methods for `LefschetzComplex` and `ConleyMorseCM`

## v0.3.10 (February 12, 2026)

- Added more methods to `sparse_show`. The first one allows two 
  additional arguments which provide labels for the rows and
  columns.
- In addition, if `cm` is a connection matrix, then the
  command `sparse_show(cm)` displays the connection matrix
  with the Conley index labels.

## v0.3.9 (September 18, 2025)

- Added extensions to `plot_planar_simplicial_morse` by Frank Pryor.

## v0.3.8 (September 13, 2025)

- Added `chain_vector` to create a sparse matrix representation of a
  chain by specifying only the cells in the support, and their 
  coefficents.
- Added `chain_support` to extract the cells in the support of a chain,
  which is given as a sparse vector.

## v0.3.7 (September 12, 2025)

- Added `forman_comb_flow` to compute Forman's combinatorial flow,
  and the associated chain homotopy.
- Added `forman_stab_flow` to compute Forman's stabilized combinatorial
  flow. Also in this case, the associated chain homotopy is returned.

## v0.3.6 (September 9, 2025)

- Added `sparse_is_identity` to check whether a sparse matrix is the 
  identity matrix.
- Added `sparse_is_equal` to test for equality of sparse matrices.
  This is also implemented as an extra method for the `==` operator.

## v0.3.5 (July 7, 2025)

This is the official archived version accompanying the JOSS submission.

## v0.3.4 (June 10, 2025)

- Added `compose_reductions`.
- Modified `lefschetz_newbasis_maps`.

## v0.3.3 (June 9, 2025)

- Added `lefschetz_cell_count` which provides cell counts in 
  each dimension.
- Added `lefschetz_newbasis` which constructs a new Lefschetz
  complex via a graded basis change.
- Added `lefschetz_newbasis_maps` which extends the previous
  function by also providing the involved isomorphisms.

## v0.3.2 (June 8, 2025)

- Added `sparse_inverse` to compute the inverse of a sparse matrix

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

