import numpy as np
from scipy.sparse import diags

def generate_safe_system(n):
    """
    Generate a linear system A x = b where A is strictly diagonally dominant,
    ensuring LU factorization without pivoting will work.

    Parameters:
        n (int): Size of the system (n x n)

    Returns:
        A (ndarray): n x n strictly diagonally dominant matrix
        b (ndarray): RHS vector
        x_true (ndarray): The true solution vector
    """

    k = [np.ones(n - 1), -2 * np.ones(n), np.ones(n - 1)]
    offset = [-1, 0, 1]
    A = diags(k, offset).toarray()

    # Solution is always all ones
    x_true = np.ones((n, 1))

    # Compute b = A @ x_true
    b = A @ x_true

    return A, b, x_true

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
    #get num rows:
    n = A.shape[0]
    
    U = A.copy()
    L = np.eye(n, dtype=np.double)
    
    for i in range(n):

    # row ops on U below i and reverse row ops for L
        factor = U[i+1:, i] / U[i, i]
        L[i+1:, i] = factor
        U[i+1:] -= factor[:, np.newaxis] * U[i]
        
    return L, U

    #arrays of zeros
    L, U = np.zeros_like(A), np.zeros_like(A)

    
    L[0, 0] = 1
    U[0, 0] = 1


#A=np.array([[4,2,0],[2,3,1],[0, 1, 2.5]])

#print(lu_factorisation(A))

def det(A):
    L, U = lu_factorisation(A)
    return np.prod(np.diag(U))

#A_large, b_large, x_large = generate_safe_system(100)
#det = det(A)
#print(f"det is {det}")


