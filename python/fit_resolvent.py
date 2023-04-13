import numpy as np
import matplotlib.pyplot as plt

# Calculates the resolvent in 1/w

nameU = "-0.10"
folder = "T0"
subfolder = ""
name_suffix = ""

file = f"data/{folder}/V_modes/{subfolder}{nameU}_resolvent{name_suffix}.txt"
one_particle = 1 / np.abs(np.loadtxt(f"data/{folder}/V_modes/{subfolder}{nameU}_one_particle{name_suffix}.txt").flatten())

M = np.loadtxt(file)
A = M[0]
B = M[1]

w_vals = 1000
w_lin = 1 / np.logspace(np.log10(1e-5), np.log10(1e-1), w_vals, dtype=complex)
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

data = -dos(w_lin).real
ax.set_xlabel(r"$\epsilon / t$")

from scipy.optimize import curve_fit
def func_ln(x, a, b):
    return a * x + b

for dat in data:
    if(dat < 0):
        print(dat)

try:
    w_log = np.log(1 / w_lin.real)
    ax.plot(w_log, np.log(data), "-", label="Real")
    popt, pcov = curve_fit(func_ln, w_log, np.log(data))
    print(f"{popt[0]} +/- {pcov[0][0]}")
    print(popt)
    ax.plot(w_log, func_ln(w_log, *popt), "--", label=r"$a \ln( x ) + b$")
    ax.set_xlabel(r"$\ln(\epsilon)$")
    ax.set_ylabel(r"$\ln(G)$")
    ax.text(-9, 7, f"$a={popt[0]}$")
    ax.text(-9, 6, f"$b={popt[1]}$")
except RuntimeError:
    print("Could not estimate curve_fit")

ax.legend()
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}_U={nameU}.pdf")
plt.show()