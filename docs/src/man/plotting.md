# Visualization

[ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl)
provides two independent visualization backends for planar complexes:
a *Luxor-based* backend that produces PDF or PNG files directly, and
a *Plots.jl-based* backend that integrates with the Julia plotting
ecosystem. The two backends are described in separate sections below.

## Luxor-Based Plotting

The functions
* [`plot_planar_simplicial`](@ref),
* [`plot_planar_simplicial_morse`](@ref),
* [`plot_planar_cubical`](@ref), and
* [`plot_planar_cubical_morse`](@ref)
produce publication-quality plots by writing directly to a PDF or PNG
file via the [Luxor.jl](https://github.com/JuliaGraphics/Luxor.jl)
package. They require an explicit output filename and a separate
coordinate vector for the vertices of the complex.

As a simple example, the following commands create a PDF image of a
small simplicial complex together with its Morse sets:

```julia
using ConleyDynamics

lc, mvf = example_forman2d()
cm = connection_matrix(lc, mvf)
coords = [[0,0],[4,0],[8,0],[2,3],[6,3],[4,6],[4,2]]
fname  = "forman2d_morse.pdf"
plot_planar_simplicial_morse(lc, coords, fname, cm.morse)
```

For cubical complexes the workflow is analogous. The function
[`get_cubical_coords`](@ref) extracts vertex coordinates from the
cube labels, so no separate coordinate vector needs to be specified
manually:

```julia
using ConleyDynamics

cc, coords = create_cubical_rectangle(4, 3)
fname = "cubical_rectangle.pdf"
plot_planar_cubical(cc, coords, fname)
```

Both families of Luxor functions also accept an
[`EuclideanComplex`](@ref) as their first argument. In that case the
coordinate vector can be omitted entirely.

## Plots-Based Plotting

The Plots.jl backend provides eight convenience functions:

| Function | Description |
|----------|-------------|
| [`plot_simplicial`](@ref) | Simplicial complex with optional Forman overlay |
| [`plot_cubical`](@ref) | Cubical complex |
| [`plot_simplicial_morse`](@ref) | Simplicial complex with colored Morse sets |
| [`plot_cubical_morse`](@ref) | Cubical complex with colored Morse sets |
| [`plot_simplicial_mvf`](@ref) | Simplicial complex with multivector field regions |
| [`plot_cubical_mvf`](@ref) | Cubical complex with multivector field regions |
| [`plot_simplicial_mv`](@ref) | Single multivector on a simplicial complex |
| [`plot_cubical_mv`](@ref) | Single multivector on a cubical complex |

### Weak Dependency

[Plots.jl](https://github.com/JuliaPlots/Plots.jl) is a *weak
dependency* of ConleyDynamics.jl. This means it is not loaded
automatically when you do `using ConleyDynamics`, and it does not
need to be installed unless you want to use the Plots.jl backend.
The eight functions above are stub definitions that become active
only after Plots.jl has been loaded. The required loading order is:

```julia
using Plots
using ConleyDynamics
```

If ConleyDynamics is loaded *before* Plots.jl, the plotting functions
will still work because Julia activates package extensions lazily.
However, loading Plots first is the recommended order.

All eight functions return a `Plots.Plot` object, so the full Plots.jl
interface applies: the result can be displayed with `display(p)`,
saved with `savefig(p, "output.png")`, or composed with other plots.

### Input: EuclideanComplex

The Plots.jl backend works exclusively with [`EuclideanComplex`](@ref)
objects, which carry embedded vertex coordinates. The simplest way to
obtain such a complex is to pass `euclidean=true` when creating a
complex:

```julia
ec, mvf = example_forman2d(euclidean=true)
cc = create_cubical_rectangle(4, 3, euclidean=true)
```

Alternatively, [`lefschetz_to_euclidean`](@ref) converts any
`LefschetzComplex` together with a coordinate vector.

### Plotting a Complex

The two basic functions display only the underlying complex, without
any dynamics:

```julia
using Plots
using ConleyDynamics

ec = create_simplicial_rectangle(5, 2, euclidean=true)
p  = plot_simplicial(ec)
display(p)
```

The keyword argument `pdim::Vector{Bool}` controls which dimensions
are drawn. The three entries correspond to vertices (dimension 0),
edges (dimension 1), and faces (dimension 2), respectively. For
example, `pdim=[false,true,true]` suppresses vertices:

```julia
p = plot_cubical(create_cubical_rectangle(4, 3, euclidean=true),
                 pdim=[false,true,true])
display(p)
```

### Plotting Morse Sets

The Morse plot functions color each Morse set in a distinct,
automatically chosen color. Cells that do not belong to any
explicitly listed Morse set are treated as implicit singletons and
rendered in additional colors.

```julia
using Plots
using ConleyDynamics

lc, mvf = example_forman2d(euclidean=true)
cm = connection_matrix(lc, mvf)
p  = plot_simplicial_morse(lc, cm.morse)
display(p)
```

When `ci=true` is passed to [`plot_simplicial_morse`](@ref), the
color of each Morse set is determined by its Conley index instead of
being chosen automatically.

For cubical complexes the usage is identical:

```julia
using Plots
using ConleyDynamics

cc, mvf = example_forman2d(euclidean=true)   # placeholder; use a cubical example
cm = connection_matrix(cc, mvf)
p  = plot_cubical_morse(cc, cm.morse)
display(p)
```

### Plotting Multivector Field Regions

The MVF plot functions render each multivector as a single colored
region. The region geometry is computed as an inflated convex hull of
the cell barycenters, so adjacent multivectors are visually separated
by a visible gap.

```julia
using Plots
using ConleyDynamics

lc, mvf = example_forman1d(euclidean=true)
p = plot_simplicial_mvf(lc, mvf)
display(p)
```

The `tubefac` parameter controls the inflation radius as a fraction
of the average edge length (default `0.05`). The `mvfalpha` parameter
sets the fill opacity (default `0.3`). Distinct colors can be assigned
to individual multivectors by passing a vector of color names:

```julia
p = plot_simplicial_mvf(lc, mvf;
                        mvfcolor=["darkorange","steelblue","seagreen"],
                        mvfalpha=0.5)
display(p)
```

By default, cells absent from the MVF are shown as implicit singleton
regions. Pass `addcritical=false` to suppress this and display only
the explicitly listed multivectors:

```julia
p = plot_simplicial_mvf(lc, mvf; addcritical=false)
display(p)
```

### Plotting a Single Multivector

The [`plot_simplicial_mv`](@ref) and [`plot_cubical_mv`](@ref)
functions highlight a single multivector on the complex background:

```julia
using Plots
using ConleyDynamics

lc, mvf = example_forman2d(euclidean=true)
mv = mvf[1]   # first multivector
p  = plot_simplicial_mv(lc, mv)
display(p)
```

Cell labels (strings) are also accepted in place of integer indices:

```julia
p = plot_simplicial_mv(lc, ["ADE", "AD", "A", "D", "E"])
display(p)
```
