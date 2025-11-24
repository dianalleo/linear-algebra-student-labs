import numpy as np
from scipy.sparse import diags
import time
import matplotlib.pyplot as plt

def system_size(A, b):

    # Validate that A is a 2D square matrix
    if A.ndim != 2:
        raise ValueError(f"Matrix A must be 2D, but got {A.ndim}D array")

    n, m = A.shape
    if n != m:
        raise ValueError(f"Matrix A must be square, but got A.shape={A.shape}")

    if b.shape[0] != n:
        raise ValueError(
            f"System shapes are not compatible: A.shape={A.shape}, "
            f"b.shape={b.shape}"
        )

    return n

def generate_safe_system(n):

    k = [np.ones(n - 1), -2 * np.ones(n), np.ones(n - 1)]
    offset = [-1, 0, 1]
    A = diags(k, offset).toarray()

    # Solution is always all ones
    x_true = np.ones((n, 1))

    # Compute b = A @ x_true
    b = A @ x_true

    return A, b, x_true

def lu_factorisation(A):

    n, m = A.shape
    if n != m:
        raise ValueError(f"Matrix A is not square {A.shape=}")

    # construct arrays of zeros
    L, U = np.zeros_like(A), np.zeros_like(A)

    # ...
    L[0, 0] = 1
    U[0, 0] = 1
    L = np.identity(n, dtype = float)

    for j in range(n):
        for i in range(j+1):
            if (i==0):
                U[i,j]= A[i,j]
                factor =U[i,i]
            elif (i<=j):
                U[i,j] = A[i,i] - L[i,j-1] * U[j-1,i]
                factor= U[i,i]

        for i in range(j+1, n):
            L[i,j]= A[i,j]/factor
    
    return L, U

def forward_substitution(A, b):
    # get size of system
    n = system_size(A, b)

    # check is lower triangular
    if not np.allclose(A, np.tril(A)):
        raise ValueError("Matrix A is not lower triangular")

    # create solution variable
    x = np.empty_like(b)

    # perform forwards solve
    for i in range(n):
        partial_sum = 0.0
        for j in range(0, i):
            partial_sum += A[i, j] * x[j]
        x[i] = 1.0 / A[i, i] * (b[i] - partial_sum)

    return x

def backward_substitution(A, b):
    # get size of system
    n = system_size(A, b)

    # check is upper triangular
    assert np.allclose(A, np.triu(A))

    # create solution variable
    x = np.empty_like(b)

    # perform backwards solve
    for i in range(n - 1, -1, -1):  # iterate over rows backwards
        partial_sum = 0.0
        for j in range(i + 1, n):
            partial_sum += A[i, j] * x[j]
        x[i] = 1.0 / A[i, i] * (b[i] - partial_sum)

    return x

def determinant(A):
    n = A.shape[0]
    L, U = lu_factorisation(A)

    det_L = 1.0
    det_U = 1.0

    for i in range(n):
        det_L *= L[i, i]
        det_U *= U[i, i]

    return det_L * det_U

def row_swap(A, b, p, q):

    # get system size
    n = system_size(A, b)
    # swap rows of A
    for j in range(n):
        A[p, j], A[q, j] = A[q, j], A[p, j]
    # swap rows of b
    b[p, 0], b[q, 0] = b[q, 0], b[p, 0]

def row_scale(A, b, p, k):

    n = system_size(A, b)

    # scale row p of A
    for j in range(n):
        A[p, j] = k * A[p, j]
    # scale row p of b
    b[p, 0] = b[p, 0] * k

def row_add(A, b, p, k, q):

    n = system_size(A, b)

    # Perform the row operation
    for j in range(n):
        A[p, j] = A[p, j] + k * A[q, j]

    # Update the corresponding value in b
    b[p, 0] = b[p, 0] + k * b[q, 0]

def gaussian_elimination(A, b, verbose=False):

    # find shape of system
    n = system_size(A, b)

    # perform forwards elimination
    for i in range(n - 1):
        # eliminate column i
        if verbose:
            print(f"eliminating column {i}")
        for j in range(i + 1, n):
            # row j
            factor = A[j, i] / A[i, i]
            if verbose:
                print(f"  row {j} |-> row {j} - {factor} * row {i}")
            row_add(A, b, j, -factor, i)

        if verbose:
            return  


A_large, b_large, x_large = generate_safe_system(100)

sizes = [2**j for j in range(1, 6)]

gausstime=[]
lutimeL=[]
lutimeU=[]

for n in sizes:
    # generate a random system of linear equations of size n
    A, b, x = generate_safe_system(n)
    L,U= lu_factorisation(A)

    # do the solve

    t0=time.time() #start time
    forward_substitution(L,b)    
    t1=time.time() #end time
    lutimeL.append(t1 - t0)

    start=time.time()
    backward_substitution(U,b)
    end=time.time()
    lutimeU.append(end - start)    

    t2=time.time() #start time
    gaussian_elimination(A,b,verbose=True)
    t3=time.time() #end time
    gausstime.append(t3 - t2)

print(determinant(A_large)) #det

plt.figure(figsize=(10, 6))
plt.plot(sizes, gausstime, marker='o', label='gaussian elimination')
plt.plot(sizes, lutimeL, marker='s', label='backward substitution')
plt.plot(sizes, lutimeU, marker='s', label='forward substitution')

plt.xlabel("matrix size")
plt.ylabel("runtime in secs")
plt.title("runtime comparison of LU factorisation vs Gaussian Elimination")
plt.show()