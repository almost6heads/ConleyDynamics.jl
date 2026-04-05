@testset "Sparse matrix creation and access" begin
    M = [1 0 1; 0 1 1; 1 1 0]
    A = sparse_from_full(M, p=2)

    @test sparse_size(A, 1) == 3
    @test sparse_size(A, 2) == 3
    @test sparse_nz_count(A) == 6

    # Round-trip
    @test full_from_sparse(A) == M

    # Entry access
    @test sparse_get_entry(A, 1, 1) == 1
    @test sparse_get_entry(A, 1, 2) == 0
    @test sparse_get_entry(A, 2, 3) == 1

    # Entry mutation
    sparse_set_entry!(A, 1, 2, 1)
    @test sparse_get_entry(A, 1, 2) == 1
    sparse_set_entry!(A, 1, 2, 0)
    @test sparse_get_entry(A, 1, 2) == 0

    # Identity and zero
    I3 = sparse_identity(3, p=2)
    @test sparse_is_identity(I3)
    @test !sparse_is_zero(I3)

    Z = sparse_zero(3, 3, p=2)
    @test sparse_is_zero(Z)
    @test !sparse_is_identity(Z)

    # Fullness/sparsity
    @test sparse_fullness(I3) + sparse_sparsity(I3) ≈ 1.0
    @test sparse_fullness(Z) == 0.0
    @test sparse_sparsity(Z) == 1.0
end

@testset "Sparse matrix arithmetic over GF(2)" begin
    A = sparse_from_full([1 0 1; 0 1 1; 1 1 0], p=2)
    B = sparse_from_full([1 1 0; 0 1 0; 0 0 1], p=2)

    # Add (= subtract over GF(2))
    C = sparse_add(A, B)
    expected = mod.([1 0 1; 0 1 1; 1 1 0] .+ [1 1 0; 0 1 0; 0 0 1], 2)
    @test full_from_sparse(C) == expected

    C2 = sparse_subtract(A, B)
    @test sparse_is_equal(C, C2)

    # Multiply by identity gives same matrix
    I3 = sparse_identity(3, p=2)
    D = sparse_multiply(I3, A)
    @test sparse_is_equal(D, A)

    # Explicit matrix multiply
    AB = sparse_multiply(A, B)
    @test full_from_sparse(AB) == mod.([1 0 1; 0 1 1; 1 1 0] * [1 1 0; 0 1 0; 0 0 1], 2)

    # Transpose
    At = sparse_transpose(A)
    @test full_from_sparse(At) == [1 0 1; 0 1 1; 1 1 0]'

    # Scale by 1 (identity operation over GF(2))
    E = sparse_scale(Int(1), A)
    @test sparse_is_equal(E, A)

    # Equality
    @test sparse_is_equal(A, A)
    @test !sparse_is_equal(A, B)
end

@testset "Sparse matrix arithmetic over rationals" begin
    Mq = [1 2; 3 4]
    Aq = sparse_from_full(Mq, p=0)

    # Round-trip
    @test full_from_sparse(Aq) == Rational{Int}.(Mq)

    # Entry type
    @test sparse_get_entry(Aq, 1, 2) == 2//1
    @test sparse_get_entry(Aq, 2, 1) == 3//1

    # Multiply by identity
    I2q = sparse_identity(2, p=0)
    @test sparse_is_equal(sparse_multiply(I2q, Aq), Aq)

    # Inverse: A * A^{-1} = I
    Aqinv = sparse_inverse(Aq)
    prod = sparse_multiply(Aq, Aqinv)
    @test sparse_is_identity(prod)

    # Scale
    Aq2 = sparse_scale(Rational{Int}(2), Aq)
    @test full_from_sparse(Aq2) == Rational{Int}.(2 .* Mq)
end

@testset "Sparse RREF and bases" begin
    # Known matrix over GF(2)
    A = sparse_from_full([1 1 0; 1 0 1; 0 1 1], p=2)
    Ar = sparse_rref(A)
    @test sparse_is_rref(Ar)

    # Non-RREF input should not be RREF
    @test !sparse_is_rref(A)

    # Kernel: 2x3 rank-2 matrix over GF(2) has 1D kernel
    G = sparse_from_full([1 0 1; 0 1 1], p=2)
    ker = sparse_basis_kernel(G)
    @test length(ker) == 1   # kernel dimension = ncols - rank = 3 - 2 = 1

    # Range: same matrix has 2D range
    rng = sparse_basis_range(G)
    @test length(rng) == 2

    # Over GF(7): example from docstring has 4-dimensional kernel
    Afull = [1 0 1 2 4 3 2; 3 4 6 2 3 4 2; 2 0 2 4 1 1 1]
    A7 = sparse_from_full(Afull, p=7)
    sparse_rref!(A7)
    ker7 = sparse_basis_kernel(A7)
    @test length(ker7) == 4
end

@testset "Sparse matrix utilities" begin
    A = sparse_from_full([1 0 1; 0 1 1; 1 1 0], p=2)

    # Permute with identity permutation gives same matrix
    perm = [1, 2, 3]
    Ap = sparse_permute(A, perm, perm)
    @test sparse_is_equal(Ap, A)

    # Permute rows and columns
    perm2 = [2, 1, 3]
    Ap2 = sparse_permute(A, perm2, perm2)
    M = full_from_sparse(A)
    @test full_from_sparse(Ap2) == M[perm2, perm2]

    # from_lists / lists_from_sparse round-trip
    nr, nc, tchar, tzero, tone, r, c, vals = lists_from_sparse(A)
    B = sparse_from_lists(nr, nc, tchar, tzero, tone, r, c, vals)
    @test sparse_is_equal(A, B)
end
