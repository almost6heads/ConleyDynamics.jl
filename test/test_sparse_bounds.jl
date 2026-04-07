@testset "Sparse matrix bounds checking" begin

    # Create a 3x4 matrix over GF(2) for testing
    sm = sparse_zero(3, 4, p=2)
    sparse_set_entry!(sm, 2, 3, 1)

    # --- sparse_get_entry / getindex ---

    # Row index out of bounds
    @test_throws BoundsError sparse_get_entry(sm, 0, 1)
    @test_throws BoundsError sparse_get_entry(sm, 4, 1)
    @test_throws BoundsError sm[0, 1]
    @test_throws BoundsError sm[4, 1]

    # Column index out of bounds
    @test_throws BoundsError sparse_get_entry(sm, 1, 0)
    @test_throws BoundsError sparse_get_entry(sm, 1, 5)
    @test_throws BoundsError sm[1, 0]
    @test_throws BoundsError sm[1, 5]

    # --- sparse_set_entry! / setindex! ---

    # Row index out of bounds
    @test_throws BoundsError sparse_set_entry!(sm, 0, 1, 1)
    @test_throws BoundsError sparse_set_entry!(sm, 4, 1, 1)
    @test_throws BoundsError (sm[0, 1] = 1)
    @test_throws BoundsError (sm[4, 1] = 1)

    # Column index out of bounds
    @test_throws BoundsError sparse_set_entry!(sm, 1, 0, 1)
    @test_throws BoundsError sparse_set_entry!(sm, 1, 5, 1)
    @test_throws BoundsError (sm[1, 0] = 1)
    @test_throws BoundsError (sm[1, 5] = 1)

    # --- sparse_get_column / sparse_get_nz_column / sparse_get_nz_row ---

    @test_throws BoundsError sparse_get_column(sm, 0)
    @test_throws BoundsError sparse_get_column(sm, 5)
    @test_throws BoundsError sparse_get_nz_column(sm, 0)
    @test_throws BoundsError sparse_get_nz_column(sm, 5)
    @test_throws BoundsError sparse_get_nz_row(sm, 0)
    @test_throws BoundsError sparse_get_nz_row(sm, 4)

    # --- sparse_add_column! ---

    sq = sparse_zero(4, 4, p=2)
    @test_throws BoundsError sparse_add_column!(sq, 0, 1, 1, 1)
    @test_throws BoundsError sparse_add_column!(sq, 1, 0, 1, 1)
    @test_throws BoundsError sparse_add_column!(sq, 5, 1, 1, 1)
    @test_throws BoundsError sparse_add_column!(sq, 1, 5, 1, 1)

    # Also test generic (non-Int) variant
    smq = sparse_zero(4, 4, p=0)  # rational field
    @test_throws BoundsError sparse_add_column!(smq, 0, 1, 1//1, 1//1)
    @test_throws BoundsError sparse_add_column!(smq, 1, 5, 1//1, 1//1)

    # --- sparse_add_row! ---

    @test_throws BoundsError sparse_add_row!(sq, 0, 1, 1, 1)
    @test_throws BoundsError sparse_add_row!(sq, 1, 0, 1, 1)
    @test_throws BoundsError sparse_add_row!(sq, 5, 1, 1, 1)
    @test_throws BoundsError sparse_add_row!(sq, 1, 5, 1, 1)

    # Also test generic variant
    @test_throws BoundsError sparse_add_row!(smq, 0, 1, 1//1, 1//1)
    @test_throws BoundsError sparse_add_row!(smq, 1, 5, 1//1, 1//1)

    # --- Valid accesses still work (regression) ---

    @test sparse_get_entry(sm, 2, 3) == 1
    @test sm[2, 3] == 1
    @test sm[1, 1] == 0

    @test sparse_get_column(sm, 3) == [0, 1, 0]
    @test sparse_get_nz_column(sm, 3) == [2]
    @test sparse_get_nz_row(sm, 2) == [3]

    sm[1, 1] = 1
    @test sm[1, 1] == 1
    sm[1, 1] = 0
    @test sm[1, 1] == 0

    # --- Internal arithmetic paths not broken (regression) ---

    A = sparse_identity(4, p=2)
    B = sparse_identity(4, p=2)
    @test sparse_is_equal(sparse_multiply(A, B), A)
    @test sparse_is_equal(sparse_add(A, A), sparse_zero(4, 4, p=2))
    @test sparse_is_equal(sparse_subtract(A, A), sparse_zero(4, 4, p=2))

    C = sparse_from_full([1 1; 0 1], p=2)
    Cinv = sparse_inverse(C)
    @test sparse_is_equal(sparse_multiply(C, Cinv), sparse_identity(2, p=2))

    D = sparse_from_full([1 2; 3 4], p=5)
    sparse_add_column!(D, 2, 1, 1, 1)
    sparse_add_row!(D, 1, 2, 1, 1)

end
