import numpy as np
import matplotlib.pyplot as plt

# Calculates the resolvent in 1/w

file = "data/resolvent.txt"

M = np.loadtxt(file)
A = M[0]
B = M[1]

w_vals = 10000
w_lin = np.linspace(-20, 20, w_vals, dtype="complex128")
w_lin += 1e-2j
off = 1

def denominator(w):
    G = 1. - w * A[len(A) - off] - w * B[len(B) - off]
    for j in range(len(A) - off - 1, -1, -1):
        G = 1. - w * A[j] - w**2 * B[j + 1] / G
    return G

def dos(w):
    G = 1. - w * A[len(A) - off] - w * B[len(B) - off]
    for j in range(len(A) - off - 1, -1, -1):
        G = 1. - w * A[j] - w**2 * B[j + 1] / G
    return -B[0] / G
    
fig, ax = plt.subplots()
#ax.plot(w_lin.real, abs(dos(w_lin).imag))
ax.plot(A[:100], 'x')
ax.plot(B[:100], 'o', mfc="none")


ax.set_xlabel(r"$\epsilon / t$")
fig.tight_layout()
plt.show()