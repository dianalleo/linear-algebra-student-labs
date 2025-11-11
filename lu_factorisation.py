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
    n, m = A.shape
    if n != m:
        raise ValueError(f"Matrix A is not square {A.shape=}")

    # construct arrays of zeros
    L, U = np.zeros_like(A), np.zeros_like(A)
    for i in range(n):
        
    #L[0,0], L[1,1], L[2,2]=1 #make L diagonal matrix

    # U[0,0]=A[0,0]
    # U[0,1]=A[0,1]
    # U[0,3]=A[0,3] #setting known U values of first row

    # L[1,0]=A[1,0]/U[0,0]
    # L[2,0]=A[2,0]

    for j in range(n):
        for i in range(j+1):
            if(i==0):
                U[i][j]=A[i][j]
            elif(i==j):
                U[i][i]=A[i][i]-L[i][j]*U[j][i]
                factor = U[i][i]
            elif(i>j):
                #
        for i in range(j+1,n):
            L[i][j]=A[i][j]/factor

    return L,U


    # ...