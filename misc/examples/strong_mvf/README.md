# Example scripts from pull request #10

These scripts were written by Justin Thorpe for pull request #10
(`ajt-math:feature/dynamic_transition_space`), and they accompany his
dissertation *Periodic Orbits in the Setting of Combinatorial Dynamics*
(George Mason University, 2026). They are kept here unchanged.

| File | Content |
|---|---|
| `example_ouroboros.jl` | The Ouroboros example: forward merging, dual sets, and dual merging on a 4x5 rectangle, with plots of every stage |
| `example_2d_circle_periodicity.jl` | A planar field with a repelling periodic orbit: multivector digraph, strongly connected components, and a split into pieces |
| `example_vdp_periodicity.jl` | The van der Pol equation on a curvature-adapted mesh |
| `example_vector_fields.jl` | The vector fields used by the scripts |

The scripts use the interface of the pull request, i.e., the functions
`create_planar_dts`, `graph_merge`, `create_dual_set`, `dual_graph_merge`,
`create_mvf`, and `digraph_propagate`, and the vector fields exported by
`example_vector_fields.jl`. These were replaced on the branch `strong-mvf`
by the following functions:

| Pull request #10 | Replacement |
|---|---|
| `create_planar_dts` | `planar_transitions` |
| `graph_merge`, `close_dts` | `create_mvf_hull` |
| `create_dual_set` | `mvf_dual_set`, `mvf_dual` |
| `create_mvf(...; strong_iso=true)`, `dual_graph_merge` | the internal cross-check `ConleyDynamics.strong_mvf_dual_merging` |
| `digraph_propagate` | `cyclic_partition_bfs` |
| `ouroboros_vf` | `example_ouroboros()` |

To run the scripts as they are, check out the merge commit of the pull
request on the branch `strong-mvf`, and include `example_vector_fields.jl`
before running the other scripts. They also need the packages `Plots`
(with the PythonPlot backend), `DifferentialEquations`, and `DistMesh`,
which are not dependencies of `ConleyDynamics.jl`.
