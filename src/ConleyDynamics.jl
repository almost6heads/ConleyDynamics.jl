"""
    module ConleyDynamics

Collection of tools for computational Conley theory.
"""
module ConleyDynamics

using Nemo
using Graphs
using Combinatorics
using Random

# Load the type definitions

include("./sparse/sparse_type.jl")
include("./cmcore/composite_types.jl")

# Load the methods

include("./cmcore/admissible_order.jl")
include("./cmcore/cm_columns.jl")
include("./cmcore/cm_create.jl")
include("./cmcore/connection_matrix.jl")
include("./cmcore/homogeneous_columns.jl")
include("./cmcore/renumber_poset.jl")
include("./cmcore/target_columns.jl")
include("./cmcore/update_low.jl")

include("./homology/ph_reduce.jl")
include("./homology/persistent_homology.jl")
include("./homology/homology.jl")
include("./homology/relative_homology.jl")

include("./utils/convert_lefschetz_sparse.jl")
include("./utils/convert_matrix_gfp.jl")
include("./utils/convert_matrix_int.jl")
include("./utils/convert_mvf.jl")
include("./utils/convert_simplices.jl")
include("./utils/create_simplicial_complex.jl")
include("./utils/permute_lefschetz_complex.jl")
include("./utils/lefschetz_closure.jl")

include("./sparse/sparse_basic_functions.jl")
include("./sparse/sparse_from_lists.jl")
include("./sparse/sparse_from_full.jl")
include("./sparse/sparse_access.jl")
include("./sparse/sparse_remove.jl")
include("./sparse/sparse_permute.jl")
include("./sparse/sparse_minor.jl")
include("./sparse/sparse_elementary_op.jl")
include("./sparse/sparse_is_sut.jl")
include("./sparse/sparse_multiply.jl")

include("./examples/example_MW_fig01.jl")
include("./examples/example_MW_fig02.jl")
include("./examples/example_MW_fig03.jl")
include("./examples/example_MW_fig04.jl")
include("./examples/example_MW_fig11.jl")
include("./examples/example_BKMW20_fig1.jl")
include("./examples/example_BKMW20_fig3.jl")
include("./examples/example_nonunique.jl")
include("./examples/example_critical_simplex.jl")
include("./examples/example_moebius.jl")

end
