# Lefschetz Complex Functions

```@meta
DocTestSetup = quote
    push!(LOAD_PATH,"../../../src/")
    using ConleyDynamics
end
```

## Lefschetz Complex Creation

```@docs
create_lefschetz_gf2
lefschetz_subcomplex
lefschetz_closed_subcomplex
lefschetz_reduction
lefschetz_reduction_maps
lefschetz_newbasis
lefschetz_newbasis_maps
compose_reductions
permute_lefschetz_complex
```

## Lefschetz Complex Queries

```@docs
lefschetz_information
lefschetz_cell_count
lefschetz_field
lefschetz_is_closed
lefschetz_is_locally_closed
validate_lefschetz_complex
```

## Topological Features

```@docs
lefschetz_boundary
lefschetz_coboundary
lefschetz_closure
lefschetz_interior
lefschetz_topboundary
lefschetz_openhull
lefschetz_lchull
lefschetz_clomo_pair
lefschetz_neighbors
lefschetz_skeleton
manifold_boundary
```

## Filters on Lefschetz Complexes

```@docs
create_random_filter
filter_shallow_pairs
filter_induced_mvf
lefschetz_filtration
lefschetz_filtration_mvf
```

## Lefschetz Helper Functions

```@docs
lefschetz_gfp_conversion
```

