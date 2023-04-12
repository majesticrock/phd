import numpy as np

M = np.loadtxt("data/test/V_modes/-0.10_gn.txt")
N = np.loadtxt("data/test/V_modes/-0.10_ng.txt")

def swap_rows(i, j):
    N[[i,j]] = N[[j,i]]

def swap_cols(i, j):
    N[:, [i,j]] = N[:, [j,i]]

basis_size = 2 * 10**2
toSwap = basis_size + np.linspace(0, basis_size, basis_size, endpoint=False, dtype=int)

for swap in toSwap:
    swap_rows(swap, swap + basis_size)

for swap in toSwap:
    swap_cols(swap, swap + basis_size)

res = np.matmul(M, N) - np.matmul(N, M)

print(np.linalg.norm(res))