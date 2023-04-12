import numpy as np

M = np.loadtxt("data/test/V_modes/-0.10_gn.txt")
N = np.loadtxt("data/test/V_modes/-0.10_ng.txt")

res = np.matmul(M, N) - np.matmul(N, M)

print(np.linalg.norm(res))