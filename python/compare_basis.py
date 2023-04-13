import numpy as np

M = np.loadtxt("data/test/V_modes/-0.10_home.txt")
N = np.loadtxt("data/test/V_modes/-0.10_work.txt")

def swap_rows(i, j):
    N[[i,j]] = N[[j,i]]

def swap_cols(i, j):
    N[:, [i,j]] = N[:, [j,i]]

basis_size = 2 * 10**2
toSwap = basis_size + np.linspace(0, basis_size, basis_size, endpoint=False, dtype=int)



res = M - N

print(np.linalg.norm(res))