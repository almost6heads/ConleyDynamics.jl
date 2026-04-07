# ConleyDynamics.jl

*Conley index and multivector fields for Julia.*

## Introduction

[ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl) is a
Julia package for studying combinatorial multivector fields using Conley
theory. The multivector fields can be studied on arbitrary Lefschetz
complexes, which include both simplicial and cubical complexes as
important special cases. The concept of combinatorial multivector field
generalizes Forman vector fields, which were originally introduced to study
Morse theory in a discrete combinatorial setting.

!!! note
    This documentation is also available in PDF format: [ConleyDynamics.pdf](ConleyDynamics.pdf).

## Features

- Data structures for Lefschetz complexes, in particular simplicial and
  cubical complexes.
- Classical Forman combinatorial vector fields and multivector fields are
  supported.
- Computation of Conley indices, connection matrices, and Conley-Morse graphs.
- Basic homology algorithms over finite fields and the rationals, including
  persistent homology and relative homology.
- Algorithms rely on a built-in sparse matrix implementation which is geared
  towards computations over finite fields and the rationals.

## Installation

To use [ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl)
please install Julia 1.10 or higher. See <https://julialang.org/downloads/>
for instructions on how to obtain Julia for your system.

At the Julia prompt simply type

```
julia> using Pkg; Pkg.add("ConleyDynamics")
```

After Julia has finished downloading and precompiling the package and
all of its dependencies, you can start using it by typing

```
julia> using ConleyDynamics
```

Whenever possible, the functions implemented in
[ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl)
have been parallelized, so as to speed up performance on multi-core
computers. This, however, is only taken advantage of if Julia is started
with multiple threads. For example, if you would like to start Julia using
8 threads, use the command

```bash
julia -t 8
```

or

```bash
julia --threads=8
```

If you want to leave the number of threads selection to the machine, 
simply use

```bash
julia --threads=auto
```

Once you are in the Julia REPL, you can see the number of used threads
using the command

```
julia> Threads.nthreads()
```

## Manual Outline

The [Tutorial](@ref) briefly explains how to get started with
[ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl).
More details, including on the underlying mathematics, are provided in
the following three sections, which cover Lefschetz complexes, homology,
and Conley theory, including connection matrices and algorithms for their
computation. After discussing a number of examples in the sections
[First Examples](@ref) and [Further Examples](@ref), the
manual concludes with a description of the sparse matrix format
underlying the package.

```@contents
Pages = [
    "man/tutorial.md",
    "man/lefschetz.md",
    "man/homology.md",
    "man/conley.md",
    "man/conley_examples.md",
    "man/examples.md",
    "man/sparse.md",
]
Depth = 2
```

## Citation

If you use
[ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl),
please cite it. There is a JOSS paper published at
[https://doi.org/10.21105/joss.08085](https://doi.org/10.21105/joss.08085),
see [wanner:25a](@cite).
The BibTeX entry for this paper is:

```bibtex
@article{wanner:25a,
   author = {Thomas Wanner},
   title = {Conley{D}ynamics.jl: {A} {J}ulia package for multivector
            dynamics on {L}efschetz complexes},
   journal = {Journal of Open Source Software},
   doi = {10.21105/joss.08085},
   volume = {10},
   number = {111},
   pages = {8085},
   url = {https://joss.theoj.org/papers/10.21105/joss.08085},
   year = {2025}
   }
```

## License 

[ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl)
is provided under an
[MIT license](https://github.com/almost6heads/ConleyDynamics.jl/blob/main/LICENSE).

