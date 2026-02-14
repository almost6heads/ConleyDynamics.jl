# Lefschetz Complex Functions

```@meta
DocTestSetup = quote
    push!(LOAD_PATH,"../../../src/")
    using ConleyDynamics
end
```

## Simplicial Complexes

```@docs
create_simplicial_complex
create_simplicial_rectangle
create_simplicial_delaunay
simplicial_torus
simplicial_klein_bottle
simplicial_projective_plane
simplicial_torsion_space
```

## Cubical Complexes

```@docs
create_cubical_complex
create_cubical_rectangle
create_cubical_box
cube_field_size
cube_information
cube_label
get_cubical_coords
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
lefschetz_skeleton
manifold_boundary
```

## Filters on Lefschetz Complexes

```@docs
create_random_filter
filter_shallow_pairs
filter_induced_mvf
```

## Lefschetz Helper Functions

```@docs
lefschetz_gfp_conversion
lefschetz_filtration
```

## Cell Subset Helper Functions

```@docs
convert_cells
convert_cellsubsets
```

## Coordinate Helper Functions

```@docs
convert_planar_coordinates
convert_spatial_coordinates
signed_distance_rectangle
signed_distance_circle
segment_intersects_rectangle
segment_intersects_circle
```

