# Plan: MVF Visualization — Stadium-Based Geometry

## Context

The current MVF plot (`plot_simplicial_mvf` / `plot_cubical_mvf`) uses inset polygons
plus bridge polygons. This has two visible problems:

1. Adjacent multivectors that share a boundary leave almost no visible gap (the inset
   amount is small and the shapes nearly touch or overlap).
2. Vertex circles have a larger radius than edge thickening, creating an inconsistent look.

The desired style (shown in `sample_mvf.png`) uses **stadium (pill) shapes** with a
single uniform tube radius `r`:
- Each vertex cell: filled circle of radius `r`
- Each edge cell: a pill-shaped tube of half-width `r` that extends from one endpoint
  (if it is in the same multivector) to the midpoint of the edge (if the other endpoint
  is in a different multivector), with a rounded cap at the "cut" end
- Each face cell: an inset polygon (vertices moved toward barycenter proportionally to `r`)
- Shapes within the same multivector overlap visually (same color = seamless union)
- The gap between adjacent multivectors appears automatically: each edge cell's tube ends
  at the midpoint with a round cap; the neighboring vertex circle is centered at the far
  endpoint; gap = edge_len/2 - 2r (positive as long as r < edge_len/4)

---

## New parameter: `tubefac`

Replace `insetfac` and `edgewidth` with a single `tubefac` parameter:

    r = tubefac * w_avg      (w_avg = average edge length across the complex)

Default: `tubefac = 0.15` gives r = 15% of average edge length.
Gap between adjacent vertex-edge multivectors = 50% - 2*15% = 20% of edge length.

---

## Files to change

| File | Change |
|------|--------|
| `src/plots/plot_data_types.jl` | Replace `insetfac::Float64` and `edgewidth::Float64` with `tubefac::Float64` in `MVFPlot`; update docstrings |
| `ext/plots_recipe_mvf.jl` | Complete rewrite of geometry helpers and recipe |

---

## Geometry helpers (new)

### `_avg_edge_length(ec)` — unchanged

### `_circle_polygon(p, r, n=16) -> (xs, ys)`
Circle of radius `r` centered at `p`, approximated as an n-gon, NaN-terminated.

### `_stadium_polygon(p1, p2, r, n_cap=8) -> (xs, ys)`
Pill (stadium) shape from `p1` to `p2` with half-width `r`.
- If `p1 == p2` (within tolerance): degenerate to a circle
- Two parallel sides + two semicircular caps (each `n_cap` polygon points)
- NaN-terminated

Construction (let `d` = unit direction from p1 to p2, `n` = left normal):

    right side:   p1 + r*n  to  p2 + r*n
    right cap at p2: arc from p2+r*n sweeping pi radians to p2-r*n  (n_cap pts)
    left side:    p2 - r*n  to  p1 - r*n
    left cap at p1: arc from p1-r*n sweeping pi radians to p1+r*n  (n_cap pts)
    close polygon, NaN

### Removed helpers
`_inset_polygon`, `_bridge_polygon`, `_shared_vertex_indices`, `_boundary_pairs_in`
— all deleted.

---

## Recipe logic (per multivector)

```julia
M_set = Set(mv)
r = data.tubefac * w_avg
col = parse(Colorant, data.mvfcolor)

for k in mv
    dim = ec.dimensions[k]

    if dim == 0
        # Circle at vertex
        xs, ys = _circle_polygon(ec.coords[k][1], r)
        @series begin
            seriestype := :shape
            fillcolor  := col
            fillalpha  --> 0.8
            linewidth  --> 0
            xs, ys
        end

    elseif dim == 1
        # Stadium: extends to v1 if v1 in M, else stops at midpoint; same for v2
        p1 = ec.coords[k][1];  p2 = ec.coords[k][2]
        mid = (p1 .+ p2) ./ 2
        bv = sparse_get_nz_column(ec.boundary, k)   # vertex cell indices
        cap1 = (length(bv) >= 1 && bv[1] in M_set) ? p1 : mid
        cap2 = (length(bv) >= 2 && bv[2] in M_set) ? p2 : mid
        xs, ys = _stadium_polygon(cap1, cap2, r)
        @series begin
            seriestype := :shape
            fillcolor  := col
            fillalpha  --> 0.8
            linewidth  --> 0
            xs, ys
        end

    elseif dim == 2
        # Inset polygon: move each vertex toward barycenter so edge gap ~= r
        cv = ec.coords[k];  nv = length(cv)
        bx = sum(p[1] for p in cv) / nv
        by = sum(p[2] for p in cv) / nv
        pts = [ begin
                    dvx = bx - p[1]; dvy = by - p[2]
                    dv  = sqrt(dvx^2 + dvy^2)
                    frac = dv > 0 ? min(r / dv, 0.45) : 0.0
                    [p[1] + frac*dvx, p[2] + frac*dvy]
                end for p in cv ]
        # polygon traversal: triangle = 1-2-3; quad = 1-3-4-2 (existing convention)
        # ... build xs, ys from pts ...
        @series begin
            seriestype := :shape
            fillcolor  := col
            fillalpha  --> 0.8
            linewidth  --> 0
            xs, ys
        end
    end
end
```

No bridge polygons are needed. The circle at a vertex and the stadium for an adjacent
edge naturally overlap (same color), merging into one connected region.

---

## Gap geometry (why it works)

For two adjacent vertex-edge multivectors {A, AB} and {B, BC} sharing vertex B
(A paired with AB, B paired with BC):

- {A, AB} emits: circle at A (radius r) + stadium from A to mid(AB) (radius r)
  - Stadium's rounded cap at mid(AB) extends to mid(AB) + r in the AB direction
- {B, BC} emits: circle at B (radius r) + stadium from B to mid(BC) (radius r)
  - Circle at B extends to B - r in the AB direction (i.e., toward A)

Gap along edge AB = |B - mid(AB)| - r - r = |AB|/2 - 2r

With r = 0.15 * |AB|:  gap = 0.5 - 0.30 = 0.20 * |AB|  (clearly visible)

---

## API changes

**`MVFPlot` struct** (`src/plots/plot_data_types.jl`):
```julia
struct MVFPlot
    complex  :: EuclideanComplex
    mvf      :: CellSubsets
    pdim     :: Vector{Bool}
    tubefac  :: Float64     # replaces insetfac + edgewidth
    mvfcolor :: String
end
```

**Convenience functions** (`ext/plots_recipe_mvf.jl`):
```julia
plot_simplicial_mvf(ec, mvf; pdim=[true,true,true], tubefac=0.15, mvfcolor="darkorange")
plot_cubical_mvf(ec,    mvf; pdim=[true,true,true], tubefac=0.15, mvfcolor="darkorange")
```

---

## Verification

```julia
using Plots, ConleyDynamics

# 1D example: vertex-edge pairs, clear pill shapes with gaps
lc, mvf = example_forman1d(euclidean=true)
p1 = plot_simplicial_mvf(lc, mvf)
display(p1)     # pills along edges, circles at vertices, ~20% edge-length gap

# 2D example
lc2, mvf2 = example_forman2d(euclidean=true)
p2 = plot_simplicial_mvf(lc2, mvf2)
display(p2)

# Larger tube
p3 = plot_simplicial_mvf(lc2, mvf2; tubefac=0.2)
display(p3)

# Custom color
p4 = plot_simplicial_mvf(lc2, mvf2; mvfcolor="mediumpurple")
display(p4)
```
