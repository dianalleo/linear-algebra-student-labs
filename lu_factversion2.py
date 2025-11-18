import numpy as np
def lu_factorisation(A):
    """
    Compute the LU factorisation of a square matrix A.

    The function decomposes a square matrix ``A`` into the product of a lower
    triangular matrix ``L`` and an upper triangular matrix ``U`` such that:

    .. math::
        A = L U

    where ``L`` has unit diagonal elements and ``U`` is upper triangular.

    Parameters
    ----------
    A : numpy.ndarray
        A 2D NumPy array of shape ``(n, n)`` representing the square matrix to
        factorise.

    Returns
    -------
    L : numpy.ndarray
        A lower triangular matrix with shape ``(n, n)`` and unit diagonal.
    U : numpy.ndarray
        An upper triangular matrix with shape ``(n, n)``.
    """
    np.asarray(A)
    n, m = A.shape
    if n != m:
        raise ValueError(f"Matrix A is not square {A.shape=}")

    # construct arrays of zeros
    L, U = np.zeros_like(A), np.zeros_like(A)

    for j in range(n):
        for i in range(j+1):
            if (i==0):
                U[i,j]= A[i,j]
                factor =U[i,i]
            elif (i<=j):
                U[i,i] = A[i,i] - L[i,j-1] * U[j-1,i]
                factor= U[i][i]

        for i in range(j+1, n):
            L[i,j]= A[i,j]/factor

    return L, U 


A=np.array([[4,2,0],[2,3,1],[0, 1, 2.5]])

print(lu_factorisation(A))