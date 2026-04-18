# Plan: Plots.jl Extension for ConleyDynamics.jl

## Context

ConleyDynamics.jl currently provides four Luxor-based plot functions that render complexes
to image files (`.pdf`, `.svg`, `.png`, `.eps`). Luxor is a full required dependency,
loaded unconditionally at module startup.

The goal is to add a parallel set of four plot functions built on Plots.jl that:
- Are activated only when the user has loaded Plots.jl before ConleyDynamics (weak dependency)
- Follow the Plots.jl recipes framework so users get full interactivity — axis manipulation,
  zooming, `savefig`, backends, etc.
- Have new names (no `_planar_` infix) so the existing Luxor functions are never touched
- Accept **only `EuclideanComplex`** — the legacy `(AbstractComplex, vertex_coords)`
  calling convention is intentionally not implemented, since it is planned for removal
  and the new functions should not carry forward that pattern
- For `EuclideanComplex`, use `ec.coords[k]` (per-cell vertex coordinates) directly —
  this handles subcomplexes correctly because each cell stores its own vertex positions
  even when the vertex cells themselves are absent from the complex
- Produce interactive display-oriented output rather than fixed-size file output

Additionally, the extension adds a **new** `plot_mvf` function family (no Luxor counterpart)
that visualizes a multivector field as shaded cell-interior regions — one connected colored
region per multivector, inset from shared boundaries so adjacent multivectors are visually
distinct.

---

## New Function Names

### Parallel to existing Luxor functions (EuclideanComplex only)

| Existing (Luxor, file output)        | New (Plots.jl, display output)   |
|--------------------------------------|----------------------------------|
| `plot_planar_simplicial`             | `plot_simplicial`                |
| `plot_planar_simplicial_morse`       | `plot_simplicial_morse`          |
| `plot_planar_cubical`                | `plot_cubical`                   |
| `plot_planar_cubical_morse`          | `plot_cubical_morse`             |

### New Plots.jl-only functions (no Luxor counterpart)

| Function                  | Purpose                                      |
|---------------------------|----------------------------------------------|
| `plot_simplicial_mvf`     | Simplicial complex with shaded MVF regions   |
| `plot_cubical_mvf`        | Cubical complex with shaded MVF regions      |

---

## Architecture

### 1. Wrapper types (main module, no Plots dependency)

Define five lightweight data-holder structs in a new file
`src/plots/plot_data_types.jl`, included from `ConleyDynamics.jl`.
These carry no dependency on Plots.jl; they exist solely to enable
`@recipe` dispatch in the extension.

Since only `EuclideanComplex` is supported, the structs hold the complex
directly. The recipes access `data.complex.coords[k]` for the per-cell
vertex coordinates of cell `k`.

```julia
# src/plots/plot_data_types.jl

export SimplicialComplexPlot, CubicalComplexPlot,
       SimplicialMorsePlot, CubicalMorsePlot,
       MVFPlot

struct SimplicialComplexPlot
    complex  :: EuclideanComplex
    mvf      :: CellSubsets       # optional Forman vector field (empty = none)
    pdim     :: Vector{Bool}      # which dimensions to draw
    labeldir :: Vector{Float64}   # vertex label directions (empty = no labels)
    labeldis :: Float64           # label offset in data-space units
end

struct CubicalComplexPlot
    complex :: EuclideanComplex
    pdim    :: Vector{Bool}
end

struct SimplicialMorsePlot
    complex   :: EuclideanComplex
    morsesets :: CellSubsets
    pdim      :: Vector{Bool}
    ci        :: Bool             # color by Conley index?
end

struct CubicalMorsePlot
    complex   :: EuclideanComplex
    morsesets :: CellSubsets
    pdim      :: Vector{Bool}
end

struct MVFPlot
    complex   :: EuclideanComplex
    mvf       :: CellSubsets
    pdim      :: Vector{Bool}
    insetfac  :: Float64          # inset fraction (default 0.2)
    edgewidth :: Float64          # half-width for edge visualization (0 = auto)
end
```

### 2. Package extension (weak dependency trigger)

Add to `Project.toml`:

```toml
[weakdeps]
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[extensions]
ConleyDynamicsPlotsExt = "Plots"

[compat]
...
Plots = "1"   # add this line
```

Create the extension module at `ext/ConleyDynamicsPlotsExt.jl`. Julia's package
extension mechanism loads it automatically when both `ConleyDynamics` and `Plots`
are loaded in the same session, and never otherwise.

```julia
# ext/ConleyDynamicsPlotsExt.jl
module ConleyDynamicsPlotsExt

using ConleyDynamics
using Plots
import RecipesBase: @recipe, @series

include("plots_recipe_simplicial.jl")
include("plots_recipe_cubical.jl")
include("plots_recipe_simplicial_morse.jl")
include("plots_recipe_cubical_morse.jl")
include("plots_recipe_mvf.jl")

end
```

### 3. Convenience functions (exported from extension)

Each extension sub-file exports one or two named functions that construct the
appropriate wrapper struct and call `plot(wrapper)`. There is exactly one method
per function — taking an `EuclideanComplex` as the first argument.

```julia
# Four parallel functions:
plot_simplicial(ec::EuclideanComplex;
                mvf::CellSubsets=CellSubsets([]),
                labeldir::Vector{<:Real}=Float64[],
                labeldis::Real=0.05,
                pdim::Vector{Bool}=[true,true,true])

plot_cubical(ec::EuclideanComplex;
             pdim::Vector{Bool}=[true,true,true])

plot_simplicial_morse(ec::EuclideanComplex,
                      morsesets::CellSubsets;
                      pdim::Vector{Bool}=[false,true,true],
                      ci::Bool=false)

plot_cubical_morse(ec::EuclideanComplex,
                   morsesets::CellSubsets;
                   pdim::Vector{Bool}=[false,true,true])

# Two new MVF functions:
plot_simplicial_mvf(ec::EuclideanComplex,
                    mvf::CellSubsets;
                    pdim::Vector{Bool}=[true,true,true],
                    insetfac::Real=0.2,
                    edgewidth::Real=0.0)

plot_cubical_mvf(ec::EuclideanComplex,
                 mvf::CellSubsets;
                 pdim::Vector{Bool}=[true,true,true],
                 insetfac::Real=0.2,
                 edgewidth::Real=0.0)
```

No `fname`, `pv`, `hfac`, `vfac`, `sfac`, or `cubefac` parameters exist in the
new API. Saving is done via `savefig(p, "file.pdf")` and display via Julia's
normal display mechanism.

---

## Recipe Design: Complex and Morse Plots

Recipes work in **data coordinates** directly (no pixel-space transformation). The
y-axis uses math orientation (y increases up). The per-cell vertex coordinates
`ec.coords[k]` are used without any transformation.

### Multi-series pattern

Each recipe uses multiple `@series` blocks — one per dimension-layer — with
NaN-separated single vectors for multi-segment paths and shapes (not one series
per cell, which would severely degrade performance for large complexes).

### SimplicialComplexPlot recipe sketch

```julia
@recipe function f(data::SimplicialComplexPlot)
    ec = data.complex

    aspect_ratio --> :equal
    legend       --> false
    framestyle   --> :none

    # Triangles (dim 2) — rendered first so edges/vertices appear on top
    if data.pdim[3]
        xs, ys = Float64[], Float64[]
        for k in ec.ncells:-1:1
            ec.dimensions[k] == 2 || continue
            cv = ec.coords[k]     # [[x1,y1],[x2,y2],[x3,y3]]
            append!(xs, [p[1] for p in cv]); push!(xs, cv[1][1], NaN)
            append!(ys, [p[2] for p in cv]); push!(ys, cv[1][2], NaN)
        end
        @series begin
            seriestype := :shape
            fillcolor  --> colorant"steelblue1"
            fillalpha  --> 0.7
            linewidth  --> 0
            xs, ys
        end
    end

    # Edges (dim 1) — NaN-separated path
    if data.pdim[2]
        xs, ys = Float64[], Float64[]
        for k in ec.ncells:-1:1
            ec.dimensions[k] == 1 || continue
            cv = ec.coords[k]     # [[x1,y1],[x2,y2]]
            push!(xs, cv[1][1], cv[2][1], NaN)
            push!(ys, cv[1][2], cv[2][2], NaN)
        end
        @series begin
            seriestype := :path
            linecolor  --> colorant"royalblue3"
            linewidth  --> 1.5
            xs, ys
        end
    end

    # Vertices (dim 0) — scatter
    if data.pdim[1]
        xs = [ec.coords[k][1][1] for k in 1:ec.ncells if ec.dimensions[k]==0]
        ys = [ec.coords[k][1][2] for k in 1:ec.ncells if ec.dimensions[k]==0]
        @series begin
            seriestype        := :scatter
            markercolor       --> colorant"royalblue4"
            markersize        --> 5
            markerstrokewidth --> 0
            xs, ys
        end
    end

    # MVF: critical cells (red scatter) + arrows (quiver series)
    # Barycenters computed as mean of ec.coords[k] entries.

    # Labels: annotations attribute with (x, y, text) tuples,
    # offset direction derived from labeldir[k].
end
```

The cubical recipe is structurally identical; quadrilateral 2-cells use four
corner points from `ec.coords[k]`.

### Arrow rendering for the `mvf` kwarg on `plot_simplicial`

Arrows use Plots.jl's `:quiver` series type. Barycenters are computed as the
mean of `ec.coords[k]` for each cell `k`. For each Forman pair `(k1, k2)` in
the MVF, an arrow is drawn from the lower-dimensional cell's barycenter to the
higher-dimensional one.

### Morse set coloring

`Colors.distinguishable_colors(n, seed; dropseed=true)` is used (Colors is already
a direct dependency). For the `ci=true` branch, `conley_index(ec, msI[m])` is
called (available via `using ConleyDynamics` inside the extension).

### Labels

Vertex labels added via the `annotations` attribute, with offsets derived from
`labeldir[k]` (direction 0–4 → E/N/W/S) scaled by `labeldis` in data-space units.

---

## Recipe Design: Multivector Field Visualization (MVFPlot)

This is **new functionality** with no Luxor counterpart. Each multivector is
rendered as a shaded geometric region — the union of the inset interiors of its
member cells plus bridge polygons across shared faces — so adjacent multivectors
are visually separated.

### Visual concept

- **Background complex**: drawn faded (low opacity) in the same style as the
  Morse plot functions
- **Each multivector**: one distinct color, filled `:shape` series
- `insetfac` (default 0.2) controls how far inward the shading reaches from each
  cell's boundary. Larger values leave a wider visible gap between multivectors.
- `edgewidth` (default auto) controls the visual thickness of the edge-cell shading.

### Inset cell geometry

For cell `k` with per-cell vertex coordinates `cv = ec.coords[k]`:

**Dimension 0 — vertex** (`cv` has 1 entry `p`):
- Small filled circle of radius `r = insetfac * w_avg / 2`, where `w_avg` is the
  average edge length in the complex (computed once per plot). Approximated as a
  regular 16-gon.

**Dimension 1 — edge** (`cv = [p1, p2]`):
- Midpoint: `m = (p1 + p2) / 2`
- Inset endpoints: `a = p1 + insetfac*(m - p1)`, `b = p2 + insetfac*(m - p2)`
- Edge direction unit vector: `d = normalize(p2 - p1)`
- Left perpendicular: `n = [-d[2], d[1]]`
- Fattened rectangle: `[a + w*n, b + w*n, b - w*n, a - w*n]`
  where `w = edgewidth` (if auto: `insetfac * |p2-p1| / 4`)

**Dimension 2 — triangle** (`cv = [p1, p2, p3]`):
- Barycenter: `b = (p1 + p2 + p3) / 3`
- Inset vertices: `pi' = pi + insetfac * (b - pi)`
- Filled triangle `[p1', p2', p3']`

**Dimension 2 — quad (cubical)** (`cv = [p1, p2, p3, p4]`):
- Barycenter: `b = (p1 + p2 + p3 + p4) / 4`
- Inset vertices: `pi' = pi + insetfac * (b - pi)`
- Filled quad `[p1', p2', p3', p4']`

### Bridge polygons (connecting adjacent cells in the same multivector)

For each boundary-incidence pair `(k_high, k_low)` within a multivector where
`k_low` is a face of `k_high`, a bridge polygon fills the gap between the two
inset regions.

**Bridge: edge `e = (v1, v2)` — triangle `t = (v1, v2, v3)`** (shared face v1-v2):

Let `n_et` be the unit normal of edge `e` pointing toward `v3`.

- Triangle-side inset corners near shared edge:
  `v1_t = v1 + insetfac*(b_t - v1)`,  `v2_t = v2 + insetfac*(b_t - v2)`
  where `b_t = (v1+v2+v3)/3`
- Edge-side rectangle corners facing the triangle (the `+w*n_et` side):
  `v1_e = a + w*n_et`,  `v2_e = b + w*n_et`
  where `a`, `b` are the inset edge endpoints
- Bridge quadrilateral: `[v1_t, v2_t, v2_e, v1_e]`

**Bridge: edge `e = (v1, v2)` — quad `q = (v1, v2, v3, v4)`** (shared face v1-v2):
- Analogous; `b_q = (v1+v2+v3+v4)/4`, same formula with 4-vertex barycenter.

**Bridge: vertex `v` — edge `e = (v, w)`**:
- The vertex circle (radius `r`) and the edge rectangle's `v`-side (at inset
  point `a`) are close enough that no explicit bridge is needed when `insetfac`
  is small. If a gap is visible, a tapered trapezoid connecting the two ±w points
  on the circle to the rectangle corners `a ± w*n` closes it.

### Assembly and rendering

For each multivector `mv`:
1. Emit one `:shape` `@series` per cell (inset polygon)
2. Emit one `:shape` `@series` per boundary-incidence pair (bridge polygon)

All series for the same multivector share the same color and alpha (~0.7). Since
overlapping `:shape` series of the same color visually merge, no polygon union
computation is required.

Color assignment: `distinguishable_colors(length(mvf_int), seed; dropseed=true)`
from Colors.jl.

### Internal helpers in `plots_recipe_mvf.jl`

```julia
_avg_edge_length(ec)                            # average edge length for auto edgewidth
_inset_polygon(ec, k, insetfac, w)              # returns (xs, ys) for cell k
_bridge_polygon(ec, k_high, k_low, insetfac, w) # returns (xs, ys) for bridge
_boundary_pairs_in(ec, mv)                      # incidence pairs within multivector mv
```

### MVFPlot recipe skeleton

```julia
@recipe function f(data::MVFPlot)
    ec        = data.complex
    mvf_int   = data.mvf isa Vector{Vector{Int}} ?
                    data.mvf : convert_cellsubsets(ec, data.mvf)
    insetfac  = data.insetfac
    w_avg     = _avg_edge_length(ec)
    edgewidth = iszero(data.edgewidth) ? insetfac * w_avg / 4 : data.edgewidth

    aspect_ratio --> :equal
    legend       --> false
    framestyle   --> :none

    # Background complex (faded, same pattern as Morse plots)
    # ...

    seed = [colorant"royalblue4", colorant"royalblue3", colorant"steelblue1"]
    cols = distinguishable_colors(length(mvf_int), seed; dropseed=true)

    for (m, mv) in enumerate(mvf_int)
        col = cols[m]
        for k in mv
            xs, ys = _inset_polygon(ec, k, insetfac, edgewidth)
            @series begin
                seriestype := :shape
                fillcolor  := col
                fillalpha  --> 0.7
                linewidth  --> 0
                xs, ys
            end
        end
        for (k_high, k_low) in _boundary_pairs_in(ec, mv)
            xs, ys = _bridge_polygon(ec, k_high, k_low, insetfac, edgewidth)
            @series begin
                seriestype := :shape
                fillcolor  := col
                fillalpha  --> 0.7
                linewidth  --> 0
                xs, ys
            end
        end
    end
end
```

---

## File Layout After Implementation

```
ConleyDynamics.jl/
├── Project.toml                             # add [weakdeps], [extensions], compat Plots
├── src/
│   ├── ConleyDynamics.jl                   # add include for plot_data_types.jl
│   └── plots/
│       ├── plot_planar_simplicial.jl        # UNCHANGED
│       ├── plot_planar_simplicial_morse.jl  # UNCHANGED
│       ├── plot_planar_cubical.jl           # UNCHANGED
│       ├── plot_planar_cubical_morse.jl     # UNCHANGED
│       └── plot_data_types.jl              # NEW: five wrapper structs
└── ext/
    ├── ConleyDynamicsPlotsExt.jl            # NEW: extension module
    ├── plots_recipe_simplicial.jl           # NEW: recipe + plot_simplicial
    ├── plots_recipe_cubical.jl              # NEW: recipe + plot_cubical
    ├── plots_recipe_simplicial_morse.jl     # NEW: recipe + plot_simplicial_morse
    ├── plots_recipe_cubical_morse.jl        # NEW: recipe + plot_cubical_morse
    └── plots_recipe_mvf.jl                 # NEW: recipe + plot_simplicial_mvf,
                                            #           plot_cubical_mvf
```

All existing files in `src/plots/` are left strictly untouched.

---

## Changes to Existing Files

| File | Change |
|------|--------|
| `Project.toml` | Add `[weakdeps]`, `[extensions]`, and `Plots = "1"` under `[compat]` |
| `src/ConleyDynamics.jl` | Add `include("./plots/plot_data_types.jl")` after type includes |
| `src/plots/plot_data_types.jl` | **New**: five wrapper structs holding `EuclideanComplex` |
| `ext/ConleyDynamicsPlotsExt.jl` | **New**: extension module, includes five sub-files |
| `ext/plots_recipe_simplicial.jl` | **New**: `@recipe` for `SimplicialComplexPlot` + `plot_simplicial` |
| `ext/plots_recipe_cubical.jl` | **New**: `@recipe` for `CubicalComplexPlot` + `plot_cubical` |
| `ext/plots_recipe_simplicial_morse.jl` | **New**: `@recipe` for `SimplicialMorsePlot` + `plot_simplicial_morse` |
| `ext/plots_recipe_cubical_morse.jl` | **New**: `@recipe` for `CubicalMorsePlot` + `plot_cubical_morse` |
| `ext/plots_recipe_mvf.jl` | **New**: `@recipe` for `MVFPlot` + `plot_simplicial_mvf`, `plot_cubical_mvf` |

---

## Verification

After implementation, verify with the following session:

```julia
using Plots
using ConleyDynamics

# --- Simplicial basic ---
sc, coords = create_simplicial_delaunay(300, 300, 30, 20)
ec = lefschetz_to_euclidean(sc, coords)
p = plot_simplicial(ec)
display(p)
savefig(p, "sc_test.pdf")

# --- Simplicial with MVF overlay (arrows for Forman pairs) ---
mvf = create_planar_mvf(sc, coords, some_vector_field)
p2 = plot_simplicial(ec; mvf=mvf)
display(p2)

# --- Subcomplex: vertex cells absent, higher-dim cells still carry coords ---
# subec = ... (some EuclideanComplex subcomplex)
# p3 = plot_simplicial(subec)   # must work correctly via ec.coords
# display(p3)

# --- Cubical ---
cc = create_cubical_rectangle(...)
ecc = lefschetz_to_euclidean(cc, ...)
p4 = plot_cubical(ecc)
display(p4)

# --- Morse sets ---
ms = morse_sets(sc, mvf)
p5 = plot_simplicial_morse(ec, ms)
display(p5)
p6 = plot_simplicial_morse(ec, ms; ci=true)
display(p6)

# --- Cubical Morse sets ---
mvf2 = create_planar_mvf(cc, coords2, some_vector_field)
ms2 = morse_sets(cc, mvf2)
p7 = plot_cubical_morse(ecc, ms2)
display(p7)

# --- MVF shaded region plot (new functionality) ---
p8 = plot_simplicial_mvf(ec, mvf)
display(p8)
p9 = plot_simplicial_mvf(ec, mvf; insetfac=0.3)
display(p9)

# --- Small example with known MVF for visual inspection ---
lc, mvf_ex = example_forman2d(euclidean=true)
p10 = plot_simplicial_mvf(lc, mvf_ex)
display(p10)

# --- Cubical MVF ---
p11 = plot_cubical_mvf(ecc, mvf2)
display(p11)

# --- Verify Luxor functions still work unchanged ---
plot_planar_simplicial(sc, coords, "lux_test.pdf")

# --- Verify extension NOT loaded without Plots (separate fresh session) ---
# using ConleyDynamics   # no Plots
# @assert !isdefined(Main, :plot_simplicial)
```

Also run the existing test suite (`] test ConleyDynamics`) to confirm no
regressions in the Luxor path.
