import numpy as np
import matplotlib.pyplot as plt
import gzip
# Calculates the resolvent in w^2

T = 0.
U = -2.0
V = -0.5

folder = "data/L=60/"
name_suffix = "SC"
type = "phase"

name = f"T={T}/U={U}_V={V}/"

file = f"{folder}{name}one_particle.dat.gz"
with gzip.open(file, 'rt') as f_open:
    one_particle = np.abs(np.loadtxt(f_open).flatten())
    roots = np.array([np.min(one_particle) * 2, np.max(one_particle) * 2])**2
    a_inf = (roots[0] + roots[1]) * 0.5
    b_inf = ((roots[1] - roots[0]) * 0.25)

file = f"{folder}{name}resolvent_{type}_{name_suffix}.dat.gz"
with gzip.open(file, 'rt') as f_open:
    M = np.loadtxt(f_open)
    A = M[0]
    B = M[1]


w_vals = 10000
w_lin = np.linspace(1e-2, 0.1, w_vals, dtype=complex)
w_lin += 1e-8j
w_lin = w_lin**2
off = 1

B_min = 1/16 * ( np.min(one_particle) - np.max(one_particle))**2 #
B_max = 1/16 * ( np.min(one_particle) + np.max(one_particle))**2 #
roots = np.array([np.min(one_particle) * 2, np.max(one_particle) * 2])

def r(w):
    ret = np.zeros(len(w), dtype=complex)
    for i in range(0, len(w)):
        root = (w[i]**2 + B_min - B_max)**2 - 4*w[i]**2 * B_min
        if(abs(w[i]) < roots[0]):
            ret[i] = (1/(w[i]*B_min)) * ( w[i]**2 + B_min - B_max + np.sqrt(root, dtype=complex) )
        elif(abs(w[i]) > roots[1]):
            ret[i] = (1/(w[i]*B_min)) * ( w[i]**2 + B_min - B_max - np.sqrt(root, dtype=complex) )
        else:
            ret[i] = (1/(w[i]*B_min)) * ( w[i]**2 + B_min - B_max - np.sqrt(root, dtype=complex) )
    return ret


def terminator(w):
    p = w - a_inf
    q = 4 * b_inf**2
    root = np.sqrt(np.real(p**2 - q), dtype=complex)
    return_arr = np.zeros(len(w), dtype=complex)
    for i in range(0, len(w)):
        if(w[i] > roots[0]):
            return_arr[i] = (p[i] - root[i]) / (2. * b_inf**2)
        else:
            return_arr[i] = (p[i] + root[i]) / (2. * b_inf**2)

    return return_arr

deviation_from_inf = np.zeros(len(A) - 1)
for i in range(0, len(A) - 1):
    deviation_from_inf[i] = abs((A[i] - a_inf) / a_inf) + abs((np.sqrt(B[i + 1]) - b_inf) / b_inf)

off_termi = len(A) - off - np.argmin(deviation_from_inf)
print("Terminating at i=", np.argmin(deviation_from_inf))
def dos(w):
    for i in range(0, len(w)):
        if(w[i].real > roots[0] and w[i].real < roots[1]):
            w[i] = w[i].real

    G = w - A[len(A) - off_termi] - B[len(B) - off_termi] * terminator( w )
    for j in range(len(A) - off_termi - 1, -1, -1):
        G = w - A[j] - B[j + 1] / G
    return B[0] / G

fig, ax = plt.subplots()

data = dos(w_lin).real
ax.set_xlabel(r"$\epsilon / t$")

from scipy.optimize import curve_fit
def func(x, a, b):
    return np.sign(x) * a * np.abs(x)**b

try:
    w_log = np.sqrt(w_lin.real)
    ax.plot(w_log, data, "-", label="Real")
    popt, pcov = curve_fit(func, w_log, data)
    print(f"{popt[0]} +/- {pcov[0][0]}")
    print(popt)
    ax.plot(w_log, func(w_log, *popt), "--", label=r"$a(x - x_0)^b$")
    ax.set_xlabel(r"$\epsilon$")
    ax.set_ylabel(r"$\Re [G]$")
    ax.text(0.45, 0.4, f"$a={popt[0]}$", transform = ax.transAxes)
    ax.text(0.45, 0.45, f"$b={popt[1]}$", transform = ax.transAxes)
    ax.set_yscale("log")
    ax.set_xscale("log")
except RuntimeError:
    print("Could not estimate curve_fit")
except ValueError:
    print("Value")

ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}_U={U}.pdf")
plt.show()