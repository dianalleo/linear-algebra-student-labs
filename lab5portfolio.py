import numpy as np

def gramSchmidt(A):
    n, m = A.shape
    Q = np.zeros((n, m))
    R = np.zeros((m, m))
    
    for i in range(m):
        v = A[:, i].copy()
        for j in range(i):
            R[j, i] = np.dot(Q[:, j], v)
            v = v - R[j, i] * Q[:, j]
        R[i, i] = np.linalg.norm(v)
        if R[i, i] > 0:
            Q[:, i] = v / R[i, i] # shld be right
    return Q, R

epsilons = [10**(-k) for k in range(6, 17)]

print(f"{'epsilon':>12} {'error1':>15} {'error2':>15} {'error3':>15}")
for e in epsilons:
    A = np.array([[1, 1 + e],
                  [1 + e, 1]], dtype=float)
    Q, R = gramSchmidt(A)

    error1 = np.linalg.norm(A - Q @ R, 2)
    error2 = np.linalg.norm(Q.T @ Q - np.eye(2), 2)
    error3 = np.linalg.norm(R - np.triu(R), 2)
    print(f"{e:12.1e} {error1:15.3e} {error2:15.3e} {error3:15.3e}")
