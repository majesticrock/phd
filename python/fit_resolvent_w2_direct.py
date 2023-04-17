import numpy as np
import matplotlib.pyplot as plt

# Calculates the resolvent in w^2

nameU = 0.1
folder = "T0"
subfolder = ""
name_suffix = "_phase_SC"

file = f"data/{folder}/V_modes/{subfolder}{nameU}_resolvent{name_suffix}.txt"

M = np.loadtxt(file)
A = M[0]
B = M[1]

w_vals = 10000
w_lin = np.linspace(0, 0.4, w_vals, dtype=complex)
w_lin = w_lin**2
w_lin += 1e-8j
off = 1

one_particle = np.abs(np.loadtxt(f"data/{folder}/V_modes/{subfolder}{nameU}_one_particle.txt").flatten())
roots = np.array([np.min(one_particle) * 2, np.max(one_particle) * 2])**2
a_inf = (roots[0] + roots[1]) * 0.5
b_inf = ((roots[1] - roots[0]) * 0.25)

def r(w):
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
off_termi = len(A) - 1 - np.argmin(deviation_from_inf)

def dos(w):
    for i in range(0, len(w)):
        if(w[i].real > roots[0] and w[i].real < roots[1]):
            w[i] = w[i].real
    G = w - A[len(A) - off_termi] - B[len(B) - off_termi] * r( w )
    for j in range(len(A) - off_termi - 1, -1, -1):
        G = w - A[j] - B[j + 1] / G
    return B[0] / G

fig, ax = plt.subplots()

data = dos(w_lin).real
ax.set_xlabel(r"$\epsilon / t$")

from scipy.optimize import curve_fit
def func(x, a, b):
    return np.sign(x-0.19865459773907004) * a * np.abs(x-0.19865459773907004)**b

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
    ax.set_yscale("symlog")
except RuntimeError:
    print("Could not estimate curve_fit")
except ValueError:
    print("Value")

ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}_U={nameU}.pdf")
plt.show()