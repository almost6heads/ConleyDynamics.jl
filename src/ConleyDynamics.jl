"""
    module ConleyDynamics

Collection of tools for computational Conley theory.
"""
module ConleyDynamics

using Graphs
using Combinatorics
using LinearAlgebra
using Random
using Printf
using Luxor
using Colors

# Load the type definitions

include("./sparse/sparse_type.jl")
include("./conley/composite_types.jl")

# Load the methods

include("./conley/cm_reduce.jl")
include("./conley/connection_matrix.jl")
include("./conley/conley_index.jl")
include("./conley/morse_sets.jl")
include("./conley/morse_interval.jl")
include("./conley/restrict_dynamics.jl")
include("./conley/remove_exit_set.jl")
include("./conley/isoinvset_information.jl")
include("./conley/forman_comb_flow.jl")
include("./conley/forman_stab_flow.jl")
include("./conley/forman_gpaths.jl")
include("./conley/chain_vector.jl")
include("./conley/chain_support.jl")

include("./homology/ph_reduce.jl")
include("./homology/persistent_homology.jl")
include("./homology/homology.jl")
include("./homology/relative_homology.jl")

include("./lefschetz/convert_cells.jl")
include("./lefschetz/convert_coordinates.jl")
include("./lefschetz/create_lefschetz_gf2.jl")
include("./lefschetz/create_simplicial_complex.jl")
include("./lefschetz/create_simplicial_rectangle.jl")
include("./lefschetz/create_simplicial_delaunay.jl")
include("./lefschetz/create_cubical_complex.jl")
include("./lefschetz/create_cubical_rectangle.jl")
include("./lefschetz/create_cubical_box.jl")
include("./lefschetz/permute_lefschetz_complex.jl")
include("./lefschetz/lefschetz_boundary.jl")
include("./lefschetz/lefschetz_openhull.jl")
include("./lefschetz/lefschetz_closure.jl")
include("./lefschetz/lefschetz_interior.jl")
include("./lefschetz/lefschetz_topboundary.jl")
include("./lefschetz/lefschetz_lchull.jl")
include("./lefschetz/lefschetz_is_closed.jl")
include("./lefschetz/lefschetz_is_locally_closed.jl")
include("./lefschetz/lefschetz_clomo_pair.jl")
include("./lefschetz/lefschetz_skeleton.jl")
include("./lefschetz/lefschetz_subcomplex.jl")
include("./lefschetz/lefschetz_closed_subcomplex.jl")
include("./lefschetz/lefschetz_filtration.jl")
include("./lefschetz/lefschetz_filtration_mvf.jl")
include("./lefschetz/lefschetz_gfp_conversion.jl")
include("./lefschetz/lefschetz_reduction.jl")
include("./lefschetz/lefschetz_reduction_maps.jl")
include("./lefschetz/lefschetz_newbasis.jl")
include("./lefschetz/lefschetz_field.jl")
include("./lefschetz/lefschetz_information.jl")
include("./lefschetz/lefschetz_cell_count.jl")
include("./lefschetz/filters.jl")
include("./lefschetz/manifold_boundary.jl")
include("./lefschetz/surfaces.jl")
include("./lefschetz/locate_planar_cellsubsets.jl")

include("./mvf/create_mvf_hull.jl")
include("./mvf/create_planar_mvf.jl")
include("./mvf/create_spatial_mvf.jl")
include("./mvf/extract_multivectors.jl")
include("./mvf/planar_nontransverse_edges.jl")
include("./mvf/mvf_information.jl")

include("./sparse/sparse_basic_functions.jl")
include("./sparse/sparse_from_lists.jl")
include("./sparse/sparse_from_full.jl")
include("./sparse/sparse_access.jl")
include("./sparse/sparse_remove.jl")
include("./sparse/sparse_permute.jl")
include("./sparse/sparse_minor.jl")
include("./sparse/sparse_elementary_op.jl")
include("./sparse/sparse_is_sut.jl")
include("./sparse/sparse_add.jl")
include("./sparse/sparse_subtract.jl")
include("./sparse/sparse_multiply.jl")
include("./sparse/sparse_scale.jl")
include("./sparse/sparse_inverse.jl")
include("./sparse/sparse_show.jl")

include("./examples/example_julia_logo.jl")
include("./examples/example_three_cm.jl")
include("./examples/example_multiflow.jl")
include("./examples/example_small_periodicity.jl")
include("./examples/example_subdivision.jl")
include("./examples/example_forman1d.jl")
include("./examples/example_forman2d.jl")
include("./examples/example_nonunique.jl")
include("./examples/example_critical_simplex.jl")
include("./examples/example_moebius.jl")
include("./examples/example_clorenz.jl")
include("./examples/example_dunce_chaos.jl")
include("./examples/example_torsion_chaos.jl")

include("./plots/plot_planar_simplicial.jl")
include("./plots/plot_planar_simplicial_morse.jl")
include("./plots/plot_planar_cubical.jl")
include("./plots/plot_planar_cubical_morse.jl")

end

