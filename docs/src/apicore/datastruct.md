# Composite Data Structures

The package relies on a number of basic composite data structures that
encompass more complicated objects. For the internal representation
of sparse matrices we refer to
[Internal Sparse Matrix Representation](@ref).

```@docs
ConleyDynamics
```

## Lefschetz Complex Type

```@docs
LefschetzComplex
Base.show(io::IO, ::MIME"text/plain", lc::LefschetzComplex)
```

## Cell Subset Types

```@docs
Cells
CellSubsets
```

## Conley-Morse Graph Type

```@docs
ConleyMorseCM
Base.show(io::IO, ::MIME"text/plain", cm::ConleyMorseCM)
```

