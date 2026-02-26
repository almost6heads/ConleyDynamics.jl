# Sparse Matrices

While Julia provides a data structure for sparse matrix computations, the
employed design decisions make it difficult to use this implementation for
computations over finite fields. This is mainly due to the fact that in
the Julia implementation, it is assumed that one can determine the zero and
one elements from the data type alone. However, a finite field data type
generally also depends on additional parameters, such as the characteristic
of the field.

Since the algorithms underlying
[ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl)
only require basic row and column operations, a specialized sparse
matrix implementation is provided in the package. It is briefly
described in the following.

## Sparse Matrix Format

Sparse matrices in this package have to be of the composite
data type [`SparseMatrix`](@ref), which is structured as follows:

```@docs; canonical=false
SparseMatrix
```

In this struct, the type `T` has to be either `Int` or `Rational{Int}`,
depending on whether the sparse matrix is interpreted as a matrix
with entries in the finite field ``GF(p)`` for some prime ``p``, or
over the field of rationals, respectively. The data type has the
following fields:

- `nrow::Int` designates the number of rows.
- `ncol::Int` gives the number of columns.
- `char::Int` specifies the characteristic of
  the underlying field ``F``. If `char=0`, then the field is the
  rationals ``\mathbb{Q}``, and one has to have `T = Rational{Int}`.
  If, on the other hand, the finite field ``F = GF(p)`` is used,
  then `char=p` has to be a prime number. In this case, the data
  type of the matrix entries has to be `T = Int`.
- `zero::T` provides 0 in the data type `T`.
- `one::T` provides 1 in the data type `T`.
- `columns::Vector{Vector{Int}}` is a vector of integer vectors, which
  contains the row indices of nonzero matrix entries in each column.
  More precisely, `columns[k]` contains an increasing list of row indices,
  which give the locations of all nonzero entries in column `k`. Note that
  the list for each colum has to be strictly increasing.
- `rows::Vector{Vector{Int}}` is a vector of integer vectors, which
  contains the column indices of nonzero matrix entries in each row.
  It is the precise dual to the previous field. This time, `rows[k]`
  contains an increasing list of column indices, which correspond to the
  nonzero entries of the matrix in the `k`-th row.
- `entries::Vector{Vector{T}}` is a vector of vectors which contains 
  the actual matrix entries. It is organized in exactly the same way as
  the field `columns`. In other words, for every `k = 1,..,ncol` the
  matrix entry in column `k` and row `columns[k][j]` is given by
  `entries[k][j]`, where `j` indexes the nonzero column entries from
  top to bottom.

This data structure is clearly redundant, in the sense that the
field `rows` is not needed to uniquely determine the matrix. However,
the type [`SparseMatrix`](@ref) is fundamental for almost every aspect
of [ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl),
as it is used to encode the incidence coefficient
map ``\kappa``, and therefore also the matrix representation of the
boundary operator ``\partial``. And for many operations on or queries
of Lefschetz complexes, one needs fast access to both the cells in
the boundary and the coboundary of a given cells. While the boundary
can easily be accessed via the field `columns`, the fast coboundary
access is aided by the field `rows`.

We would like to point out that in view of the different underlying
fields, sparse matrices should only be manipulated using the specific
commands provided by the package. These are described in detail below.
If there is a need for additional functionality beyond these first
methods, it can be added at a later point in time.

## Creating Sparse Matrices

The package provides a number of methods for creating sparse 
matrices with the data type [`SparseMatrix`](@ref). These are
geared towards their usage within
[ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl)
and are therefore by no means exhaustive:

- [`sparse_from_full`](@ref) is usually invoked in the form
  `A = sparse_from_full(AF, p=PP)`. The first input argument `AF`
  has to be a regular Julia integer matrix. This matrix is then
  converted to sparse format and returned as `A`. If the optional
  parameter `p` is omitted, the resulting sparse matrix is over
  the rational numbers ``\mathbb{Q}``, otherwise it is over the
  finite field with characteristic `PP`.
- [`full_from_sparse`](@ref) converts a given sparse matrix into
  a standard full matrix in Julia. The data type of the entries is
  either `Rational{Int}` or `Int`, depending on whether the sparse
  input matrix is considered over the rationals ``\mathbb{Q}`` or
  over a finite field, respectively. When invoking this command, 
  be mindful of the size of the sparse matrix!
- [`sparse_from_lists`](@ref) creates a sparse matrix solely
  based on its nonzero entries and their locations. It expects
  the following required input arguments, in the order they are
  listed:
  - `nr::Int`: Number of rows
  - `nc::Int`: Number of columns
  - `tchar`: Field characteristic, which has to be 0 if
    ``F = \mathbb{Q}`` and a positive prime otherwise
  - `tzero::T`: Number 0 of type `T`
  - `tone::T`:  Number 1 of type `T`
  - `r::Vector{Int}`: Vector of row indices
  - `c::Vector{Int}`: Vector of column indices
  - `v::Vector{T}`: Vector of matrix entries
  The function assumes that the vectors `r`, `c`, and `v` have
  the same length and that the matrix has entry `v[k]` at the
  location `(r[k],c[k])`. Zero entries will be ignored, and multiple
  entries for the same matrix position raise an error. Furthermore,
  if `tchar>0`, then the entries in `v` are all replaced by their
  values modulo `tchar`. As mentioned before, if `tchar=0` then
  the entry type has to be `T = Rational{Int}`, otherwise
  we have `T = Int`.
- [`lists_from_sparse`](@ref) takes a sparse matrix and
  disassembles it into the separate ingredients specified
  in the discussion of the previous function. In this sense, it
  is precisely the inverse method of [`sparse_from_lists`](@ref).
- [`sparse_identity`](@ref) creates a sparse identity matrix.
  It is invoked as `A = sparse_identity(n, p=PP)`, and returns
  a sparse identity matrix `A` with `n` rows and `n` columns.
  If the optional characteristic parameter specified and positive,
  then the matrix is considered over the finite field with
  characteristic `PP`, otherwise it is over the rationals
  ``\mathbb{Q}``.
- [`sparse_zero`](@ref) creates a sparse zero matrix.
  It is invoked as `A = sparse_zero(nr, nc, p=PP)`, and returns
  a sparse zero matrix `A` with `nr` rows and `nc` columns.
  If the optional characteristic parameter specified and positive,
  then the matrix is considered over the finite field with
  characteristic `PP`, otherwise it is over the rationals
  ``\mathbb{Q}``.

Of these methods, the function [`sparse_from_lists`](@ref)
provides the easiest and quickest way to create a sparse
matrix.

## Sparse Matrix Access

Access to the entries of sparse matrices is provided via the
following commands:

- [`sparse_get_entry`](@ref) extracts the matrix entry `val`
  of the matrix `A` located in row `ri` and column `ci`, if it
  is invoked using the command `val = sparse_get_entry(A,ri,ci)`.
- [`sparse_set_entry!`](@ref) sets the matrix entry of the matrix
  `A` located in row `ri` and column `ci` to the value 'val', if it
  is invoked using the command `sparse_set_entry!(A,ri,ci,val)`.
  Internally, this commands makes sure that the above-defined
  format of the fields of a sparse matrix is preserved. Note that
  the data type of `val` has to match the type of `A.zero`. If
  the matrix is a sparse matrix over the integers, it is considered
  as a sparse matrix over a finite field. Thus, the value `val`
  is automatically reduced modulo `A.char` before being assigned
  to the entry.
- [`sparse_get_column`](@ref) is invoked as
  `Acol = sparse_get_column(A,ci)`, and it returns the full `ci`-th
  column of the matrix `A` as a `Vector{T}` of length `A.nrow`.
- [`sparse_get_nz_column`](@ref) returns the row indices for the
  nonzero entries in the `ci`-th column of the sparse matrix `A`,
  if invoked as `rivec = sparse_get_nz_column(A,ci)`.
- [`sparse_get_nz_row`](@ref) returns the column indices for the
  nonzero entries in the `ri`-th row of the sparse matrix `A`,
  if invoked as `civec = sparse_get_nz_row(A,ri)`.
- [`sparse_minor`](@ref) creates a minor from a given sparse
  matrix `A`. For this, one needs to specify the row and column
  indices of the minor in the integer vectors `rvec` and `cvec`,
  respectively, and then invoke the function using the command
  `AM = sparse_minor(A,rvec,cvec)`. Note that the entries in
  `rvec` and `cvec` do not have to be in increasing order, but
  they are not allowed to contain repeated indices.

One can also read and set sparse matrix values using the overloaded
methods `y = A[i,j]` and `A[i,j] = val`. In the latter case, if the
matrix is defined as `SparseMatrix{Int}`, then it is interpreted
as being over a finite field. In this case, the value `val` is 
automatically reduced via modular arithmetic modulo `A.char` before
the assignment.

## Elementary Matrix Operations

The following commands perform the basic sparse matrix operations
that are needed for the functionality of the package:

- [`sparse_add_column!`](@ref) is invoked using the form
  `sparse_add_column!(A,ci1,ci2,cn,cd)`, and it replaces the
  `ci1`-th column `column[ci1]` of `A` by `column[ci1] +
  (cn/cd) * column[ci2]`. This operation automatically performs
  the computations over the field ``F`` underlying the sparse
  matrix `A`. In other words, if this field is finite, then it
  determines the inverse of the argument `cd` as part of the
  computation.
- [`sparse_add_row!`](@ref) is invoked using the form
  `sparse_add_row!(A,ri1,ri2,cn,cd)`, and it replaces the
  `ri1`-th row `row[ri1]` of `A` by `row[ri1] + (cn/cd) *
  row[ri2]`. As before, this operation automatically performs
  the computations over the field ``F`` underlying the sparse
  matrix `A`.
- [`sparse_permute`](@ref) creates a new sparse matrix by
  permuting the row and column indices. It is invoked using the
  command `AP = sparse_permute(A,pr,pc)`, and the integer
  vectors `pr` and `pc` have to describe the row and column
  permutations, respectively.
- [`sparse_remove!`](@ref) is invoked as `sparse_remove!(A,ri,ci)`
  and removes the sparse matrix entry in the `ri`-th row
  and `ci`-th colum, i.e., it effectively sets the entry
  equal to zero.
- [`sparse_add`](@ref) computes the matrix sum of 
  two sparse matrices. Exceptions are raised if the matrix
  sum is not defined, or if the involved sparse matrices
  are defined over different fields. One can also use the
  operator form `A+B` to compute the sum of sparse matrices.
- [`sparse_subtract`](@ref) computes the matrix difference of 
  two sparse matrices. Exceptions are raised if the matrix
  difference is not defined, or if the involved sparse matrices
  are defined over different fields. One can also use the
  operator form `A-B` to compute the sum of sparse matrices.
- [`sparse_multiply`](@ref) computes the matrix product of 
  two sparse matrices. Exceptions are raised if the matrix
  product is not defined, or if the involved sparse matrices
  are defined over different fields. One can also use the
  operator form `A*B` to compute the product of sparse matrices.
- [`sparse_scale`](@ref) computes the scalar product of a
  number and a sparse matrix. An exception is raised if the
  scalar and the matrix entries are not of the same type.
  One can also use the operator form `sfac*A` to compute the
  scalar product.
- [`sparse_inverse`](@ref) computes the inverse matrix of 
  a sparse matrix. Exceptions are raised if the matrix
  is not square or not invertible.

There are also the useful helper functions [`scalar_inverse`](@ref),
[`scalar_multiply`](@ref), and [`scalar_add`](@ref) which compute
the inverse of a scalar, the product of two scalars, or the sum of
two scalars, respectively. Depending on the input type, these
computations are performed either over the rationals, or in
modular arithmetic. See the function documentations for more
details.

As mentioned earlier, additional operations can easily be
implemented if they become necessary.

## Sparse Matrix Information

Finally,
[ConleyDynamics.jl](https://almost6heads.github.io/ConleyDynamics.jl)
provides the following functions 
for quickly extracting certain information from sparse matrices:

- [`sparse_size`](@ref) is invoked as `size = sparse_size(A,dim)`,
  and it returns the number of rows if `dim=1`, or the number of
  columns for `dim=2`.
- [`sparse_low`](@ref) returns the largest row index `ri` of a
  nonzero entry in the `ci`-th column of the matrix `A`, if used
  in the form `ri = sparse_low(A,ci)`. In other words, it returns
  the row index of the lowest nonzero matrix entry in the column.
- [`sparse_is_zero`](@ref) checks whether a sparse matrix is the
  zero matrix.
- [`sparse_is_identity`](@ref) checks whether a sparse matrix is the
  identity matrix.
- [`sparse_is_equal`](@ref) checks whether two sparse matrices are
  the same. For this, they not only have to have the same size
  and the same entries, they also need to be defined over the same
  field. This function can also be invoked using `A == B`.
- [`sparse_is_sut`](@ref) checks whether a given sparse matrix is
  strictly upper triangular, and returns either `true` or `false`.
- [`sparse_fullness`](@ref) returns the fullness of a sparse
  matrix as a floating point number. Here fullness refers to
  the ratio of the number of nonzero matrix elements and the
  total number of matrix entries.
- [`sparse_sparsity`](@ref) computes the sparseness of a sparse
  matrix, which is defined as ``1`` minus its fullness, i.e., it
  is the ratio of the number of zero matrix elements and the total
  number of matrix entries.
- [`sparse_nz_count`](@ref) determines the number of nonzero
  matrix entries.
- [`sparse_show`](@ref) can be used to display a sparse matrix
  in traditional matrix form at the Julia REPL prompt. This function
  has additional methods which allow the user to display 
  row and column labels, which have to be specified as second and third
  arguments. In addition, if `cm` is a connection matrix, then the
  command `sparse_show(cm)` shows the connection matrix including the
  Conley index labels.

