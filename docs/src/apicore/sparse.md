# Sparse Matrix Functions

## Internal Sparse Matrix Representation

```@docs
SparseMatrix
```

## Access Functions

```@docs
sparse_get_entry
Base.getindex(matrix::SparseMatrix, ri::Int, ci::Int)
sparse_set_entry!
Base.setindex!(matrix::SparseMatrix, val, ri::Int, ci::Int)
sparse_get_column
sparse_get_nz_column
sparse_get_nz_row
sparse_minor
```

## Basic Functions

```@docs
sparse_size
sparse_low
sparse_is_zero
sparse_is_identity
sparse_is_equal
Base.:(==)(::SparseMatrix,::SparseMatrix)
sparse_is_rref
sparse_is_sut
sparse_fullness
sparse_sparsity
sparse_nz_count
sparse_show
Base.show(io::IO, ::MIME"text/plain", sm::SparseMatrix)
```

## Conversion Functions

```@docs
sparse_from_full
full_from_sparse
sparse_from_lists
lists_from_sparse
```

## Matrix Creation

```@docs
sparse_identity
sparse_zero
sparse_diagonal
sparse_hcat
sparse_vcat
sparse_hvcat
sparse_hvncat
sparse_cat
Base.hcat(::SparseMatrix...)
Base.vcat(::SparseMatrix...)
Base.hvcat(::Tuple{Vararg{Int}}, ::SparseMatrix...)
Base.hvcat(::Int, ::SparseMatrix...)
Base.hvncat(::Tuple{Int,Int}, ::Bool, ::SparseMatrix...)
```

## Elementary Matrix Operations

```@docs
sparse_add_column!
sparse_add_row!
sparse_rref!
sparse_rref
sparse_solve
sparse_basis_kernel
sparse_basis_range
sparse_permute
sparse_transpose
Base.adjoint(::SparseMatrix)
sparse_inverse
sparse_remove!
sparse_add
sparse_subtract
sparse_multiply
sparse_scale
Base.:+(::SparseMatrix,::SparseMatrix)
Base.:-(::SparseMatrix,::SparseMatrix)
Base.:*(::SparseMatrix,::SparseMatrix)
Base.:*(::Int,::SparseMatrix)
Base.:*(::Rational{Int},::SparseMatrix)
Base.:*(::Any,::SparseMatrix)
```

## Sparse Helper Functions

```@docs
scalar_inverse
scalar_multiply
scalar_add
```

