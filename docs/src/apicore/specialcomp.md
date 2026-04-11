# Special Complex Functions

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

## Cell Subset Helper Functions

```@docs
convert_cells
convert_cellsubsets
cellsubsets_to_cells
cellsubset_bounding_box
cellsubset_distance
cellsubset_planar_area
locate_planar_cellsubsets
cellsubset_location_rectangle
cellsubset_location_circle
```

## Geometry Helper Functions

```@docs
convert_planar_coordinates
convert_spatial_coordinates
signed_distance_rectangle
signed_distance_circle
segment_intersects_rectangle
segment_intersects_circle
```

