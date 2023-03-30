import numpy as np
import matplotlib.pyplot as plt

# Calculates the resolvent in 1/w

nameU = "-0.10"
folder = "T0"
subfolder = ""

file = f"data/{folder}/V_modes/{subfolder}{nameU}_resolvent.txt"
one_particle = 1 / np.abs(np.loadtxt(f"data/{folder}/V_modes/{subfolder}{nameU}_one_particle.txt").flatten())

M = np.loadtxt(file)
A = M[0]
B = M[1]

w_vals = 10000
w_lin = 1 / np.linspace(1e-6, 1e-1, w_vals, dtype=complex)
w_lin += 1e-8j
off = 1

B_min = 1/16 * ( np.min(one_particle) - np.max(one_particle))**2 #
B_max = 1/16 * ( np.min(one_particle) + np.max(one_particle))**2 #
roots = np.array([np.sqrt((np.sqrt(B_min) - np.sqrt(B_max))**2), np.sqrt((np.sqrt(B_min) + np.sqrt(B_max))**2)])

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

def r2(w):
    ret = w - A[-2]
    for j in range(0, 1000):
        for k in range(len(A) - 1, -1, len(A) - 4):
            ret = w - A[k] - B[k] / ret

    return ret

def dos_r(w):
    G = w - A[len(A) - off] - B[len(B) - off] * r( w )
    for j in range(len(A) - off - 1, -1, -1):
        G = w - A[j] - B[j + 1] / G
    return w * B[0] / G

def dos(w):
    G = w - A[len(A) - off] - B[len(B) - off] #* r( w )
    for j in range(len(A) - off - 1, -1, -1):
        G = w - A[j] - B[j + 1] / G
    return w * B[0] / G
    
fig, ax = plt.subplots()
#ax.plot(1 / w_lin.real, -dos( w_lin ).imag, "-", label="Lanczos 200")
#ax.plot(1 / w_lin.real, -dos( w_lin ).imag, "--", label="Lanczos 200, Ter")
#R = np.loadtxt(f"data/{folder}/V_modes/{subfolder}{nameU}.txt")
#ax.plot(np.linspace(-10, 10, len(R)), R, "--", label="Exact")
data = -dos(w_lin).real
ax.plot(1 / w_lin.real, data, "-", label="Real")

ax.set_yscale("log")
ax.set_xscale("log")

from scipy.optimize import curve_fit
def func(x, a, b, c):
    return a * (x-c)**(-b)
popt, pcov = curve_fit(func, 1 / w_lin.real, data, (0.165, 2, 1e-7))
print(f"{popt[1]} +/- {pcov[1][1]}")
print(popt)
ax.plot(1 / w_lin.real, func(1 / w_lin.real, *popt), "--", label="$a ( x - c )^{-b}$")
ax.text(3e-6, 1e5, f"$a={popt[0]}$")
ax.text(3e-6, 1e4, f"$b={popt[1]}$")
ax.text(3e-6, 1e3, f"$c={popt[2]}$")

#ax.plot(A, 'x', label="$a_i$")
#ax.plot(B, 'o', label="$b_i$")
#ax.set_yscale("symlog")
#ax.set_ylim(-10, 10)

ax.legend()
ax.set_xlabel(r"$\epsilon / t$")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}_U={nameU}.pdf")
plt.show()