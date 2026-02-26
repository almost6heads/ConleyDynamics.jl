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
sparse_is_sut
sparse_identity
sparse_zero
sparse_fullness
sparse_sparsity
sparse_nz_count
sparse_show
Base.show(io::IO, ::MIME"text/plain", sm::SparseMatrix)
```

## Elementary Matrix Operations

```@docs
sparse_add_column!
sparse_add_row!
sparse_permute
sparse_inverse
sparse_remove!
sparse_add
sparse_subtract
sparse_multiply
sparse_scale
Base.:+(::SparseMatrix,::SparseMatrix)
Base.:-(::SparseMatrix,::SparseMatrix)
Base.:*(::SparseMatrix,::SparseMatrix)
Base.:*(::Any,::SparseMatrix)
```

## Conversion Functions

```@docs
sparse_from_full
full_from_sparse
sparse_from_lists
lists_from_sparse
```

## Sparse Helper Functions

```@docs
scalar_inverse
scalar_multiply
scalar_add
```

