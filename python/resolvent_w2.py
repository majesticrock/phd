import numpy as np
import matplotlib.pyplot as plt

# Calculates the resolvent in w^2

nameU = "-0.10"
folder = "test"
subfolder = ""
name_suffix = ""

file = f"data/{folder}/V_modes/{subfolder}{nameU}_resolvent{name_suffix}.txt"
one_particle = np.abs(np.loadtxt(f"data/{folder}/V_modes/{subfolder}{nameU}_one_particle{name_suffix}.txt").flatten())

M = np.loadtxt(file)
A = M[0]
B = M[1]

w_vals = 10000
w_lin = np.linspace(0, 10, w_vals, dtype=complex)**2
w_lin += 1e-2j
#w_lin = w_lin**2
off = 1

roots = np.array([np.min(one_particle) * 2, np.max(one_particle) * 2])**2
a_inf = (roots[0] + roots[1]) * 0.5
b_inf = ((roots[0] - roots[1]) * 0.25)

def r(w):
    w = np.asarray(w)
    scalar_input = False
    if w.ndim == 0:
        w = w[None]  # Makes w 1D
        scalar_input = True

    p = w - a_inf
    q = 4 * b_inf**2
    root = np.sqrt(np.real(p**2 - q), dtype=complex)
    return_arr = np.zeros(len(w), dtype=complex)
    for i in range(0, len(w)):
        if(w[i] > roots[0]):
            return_arr[i] = (p[i] - root[i]) / (2. * b_inf**2)
        else:
            return_arr[i] = (p[i] + root[i]) / (2. * b_inf**2)

    if scalar_input:
        return np.squeeze(return_arr)
    return return_arr

off_termi = 100
def dos_r(w):
    G = w - A[len(A) - off_termi] - B[len(B) - off_termi] * r( w )
    for j in range(len(A) - off_termi - 1, -1, -1):
        G = w - A[j] - B[j + 1] / G
    return B[0] / G

def dos(w):
    G = w - A[len(A) - off] - B[len(B) - off] #* r( w )
    for j in range(len(A) - off - 1, -1, -1):
        G = w - A[j] - B[j + 1] / G
    return B[0] / G
    
fig, ax = plt.subplots()
ax.plot(np.sqrt(w_lin.real), -dos( w_lin ).imag, "-", label="Lanczos 200")
ax.plot(np.sqrt(w_lin.real), -dos_r( w_lin ).imag, "--", label="Lanczos Termi")
ax.axvspan(np.sqrt(roots[1]), np.sqrt(roots[0]), alpha=.2, color="purple", label="Continuum")
#ax.plot(np.sqrt(w_lin.real), -r(w_lin).imag, "--", label="Im T(w)")
#ax.plot(np.sqrt(w_lin.real), r(w_lin).real, "--", label="Re T(w)")
#R = np.loadtxt(f"data/{folder}/V_modes/{subfolder}{nameU}{name_suffix}.txt")
#ax.plot(np.linspace(0, 10, len(R)), R, "--", label="Exact")
#ax.plot(np.sqrt(w_lin.real), dos(w_lin).real, ":", label="Real")

#print(a_inf, abs(b_inf))
#ax.axhline(a_inf, color="k", label="$a_\\infty$")
#ax.axhline(abs(b_inf), color="k", linestyle="--", label="$b_\\infty$")
#ax.plot(A, 'x', label="$a_i$")
#ax.plot(np.sqrt(B), 'o', label="$b_i$")
#ax.set_yscale("log")
#ax.set_xscale("log")
#ax.set_ylim(-10, 10)

ax.legend()
ax.set_xlabel(r"$\epsilon / t$")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}_U={nameU}.pdf")
plt.show()
