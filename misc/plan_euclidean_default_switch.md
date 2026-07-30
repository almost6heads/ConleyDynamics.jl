# Plan: Switching the `euclidean` keyword default to `true`

Status: on hold, audit only — nothing implemented yet. Recorded here so
the investigation doesn't need to be redone.

## Motivation

Several Lefschetz-complex-creation functions accept an `euclidean::Bool`
keyword controlling whether the result carries embedded coordinates
(`EuclideanComplex`) or not (`LefschetzComplex` + a separate coords
value). All of them currently default to `false`. Considering flipping
the default to `true` throughout.

## Audit: every function with an `euclidean` keyword

All default to `false` today. Two families, with materially different
risk profiles (see below).

| Function | File | Default | Family |
|---|---|---|---|
| `create_simplicial_complex` | `lefschetz/create_simplicial_complex.jl` | *(no `euclidean` kwarg — always returns `EuclideanComplex`)* | n/a |
| `create_simplicial_rectangle` | `lefschetz/create_simplicial_rectangle.jl` | `false` | `create_*` (arity changes) |
| `create_simplicial_delaunay` (both methods) | `lefschetz/create_simplicial_delaunay.jl` | `false` | `create_*` (arity changes) |
| `create_cubical_complex` | `lefschetz/create_cubical_complex.jl` | `false` | `create_*` (arity changes) |
| `create_cubical_rectangle` | `lefschetz/create_cubical_rectangle.jl` | `false` | `create_*` (arity changes) |
| `create_cubical_box` | `lefschetz/create_cubical_box.jl` | `false` | `create_*` (arity changes) |
| `example_forman1d` | `examples/example_forman1d.jl` | `false` | `example_*` (type only) |
| `example_forman2d` | `examples/example_forman2d.jl` | `false` | `example_*` (type only) |
| `example_nonunique` | `examples/example_nonunique.jl` | `false` | `example_*` (type only) |
| `example_three_cm` | `examples/example_three_cm.jl` | `false` | `example_*` (type only) |

## Finding 1: `Base.show` is a non-issue

`src/conley/composite_types.jl` already has parallel `Base.show(::IO,
::MIME"text/plain", ...)` methods for both `LefschetzComplex` and
`EuclideanComplex` — same layout (type header, field list, dim, ncells,
boundary line), differing only in the type name (correctly, via
`typeof`) and `EuclideanComplex` listing the extra `coords` field.
No doctest anywhere captures this literal text, and the only
`typeof(...)` doctest checks (`example_forman1d.jl`, `example_forman2d.jl`)
apply to the already-explicit `euclidean=true` call, not the default
path. Nothing to change here.

## Finding 2: the `example_*` family is low-risk

`example_forman1d`, `example_forman2d`, `example_nonunique`,
`example_three_cm` return the *same tuple arity* regardless of
`euclidean` — only the type of the complex element changes
(`LefschetzComplex` → `EuclideanComplex`). No test or doctest asserts
on that type for the default (unqualified) call path. Flipping these
four defaults independently should be safe.

## Finding 3 (the blocker): the `create_*` family changes return arity, not just type

This is more serious than a display/doctest concern — it's a hard
runtime crash. The five `create_*` functions (`create_simplicial_rectangle`,
`create_cubical_rectangle`, `create_cubical_box`, `create_cubical_complex`,
both `create_simplicial_delaunay` methods) all follow this pattern:

```julia
if euclidean
    return lefschetz_to_euclidean(lc, coords)   # 1 value: EuclideanComplex
else
    return lc, coords                            # 2 values: (LefschetzComplex, coords)
end
```

Every existing call site destructures the default (`euclidean=false`)
path as `x, coords = create_XXX(...)`. Flipping the default breaks all
of them with:

```
MethodError: no method matching iterate(::EuclideanComplex)
```

(reproduced directly against `create_simplicial_rectangle(5,2,
euclidean=true)` on 2026-07-30). This is not cosmetic — it's a crash,
not a wrong-output test failure.

Known call sites that would break (found via grep, not necessarily
exhaustive):
- `test/runtests.jl`: `create_simplicial_rectangle(5,2)`,
  `create_cubical_complex(cubes)`, `create_cubical_rectangle(5,2,
  randomize=0.2)`.
- `src/mvf/create_planar_mvf.jl`, `src/mvf/create_spatial_mvf.jl`.
- `src/plots/plot_planar_simplicial.jl`, `src/plots/plot_planar_cubical.jl`.
- Numerous unqualified calls across `docs/src/man/*.md`
  (`tutorial.md`, `lefschetz.md`, `homology.md`, `flowexamples.md`,
  `plotting.md`).

## Options for the `create_*` family (not yet decided)

1. **Make the branches arity-consistent first**, e.g. always return
   `(complex, coords)` where `coords` is either the raw vector (Lefschetz
   case) or extracted from `complex.coords` (Euclidean case) — then the
   default flip becomes purely a type change, same low risk as the
   `example_*` family. Changes the `euclidean=true` call signature too
   (currently single-value), so still a breaking change for *existing*
   `euclidean=true` callers, but a much smaller, mechanical one.
2. **Leave arity as-is, flip the default anyway**, and fix every call
   site (library + 15+ manual call sites) to match. Straightforward but
   touches a lot of files, including prose in the manual chapters that
   explains the two-tuple return.
3. **Don't flip `create_*` defaults**, only flip the `example_*` ones
   (Finding 2), since those are safe in isolation and still move in the
   direction the user wants for the example gallery.

## Next steps when resumed

- Decide between the three options above (leaning toward option 1 if
  the goal is really "always embed by default" everywhere, since it
  also cleans up the arity inconsistency independent of the default
  value).
- Once decided, sweep and update all call sites listed above (library
  code, `test/runtests.jl`, all `docs/src/man/*.md` examples using these
  functions, and their docstrings' "When `euclidean=false`
  (default).../When `euclidean=true`..." prose).
- Re-run `Pkg.test()` and a full `make.jl --local-html` doctest pass
  after any change, same as was done for the `ap/` migration and its
  API docs.
