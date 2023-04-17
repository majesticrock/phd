import numpy as np
from scipy.optimize import *
from scipy.misc import derivative
import matplotlib.pyplot as plt

names = np.array([0.005, 0.01, 0.05, 0.1, 0.2, 0.35, 0.5])
folder = "T0"
subfolder = ""
name_suffix = "SC"
fig, ax = plt.subplots()
lss = ["-", "--"]
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

peak_weights = np.zeros(len(names))
peak_positions = np.zeros(len(names))

for jdx, nameU in enumerate(names):
    one_particle = np.abs(np.loadtxt(f"data/{folder}/V_modes/{subfolder}{nameU}_one_particle.txt").flatten())
    roots = np.array([np.min(one_particle) * 2, np.max(one_particle) * 2])**2
    a_inf = (roots[0] + roots[1]) * 0.5
    b_inf = ((roots[1] - roots[0]) * 0.25)
    
    #ax.axvspan(np.sqrt(roots[1]), np.sqrt(roots[0]), alpha=.2, color="purple", label="Continuum")
    file = f"data/{folder}/V_modes/{subfolder}{nameU}_resolvent_higgs_{name_suffix}.txt"

    M = np.loadtxt(file)
    A = M[0]
    B = M[1]

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

    deviation_from_inf = np.zeros(len(A) - 1)
    for i in range(0, len(A) - 1):
        deviation_from_inf[i] = abs((A[i] - a_inf) / a_inf) + abs((np.sqrt(B[i + 1]) - b_inf) / b_inf)

    off_termi = len(A) - 1 - np.argmin(deviation_from_inf)
    def denominator(w):
        w = w**2
        G = w - A[len(A) - off_termi] - B[len(B) - off_termi] * r( w )
        for j in range(len(A) - off_termi - 1, -1, -1):
            G = w - A[j] - B[j + 1] / G
        return G / B[0]
    
    def abs_den(w):
        return np.abs(denominator(w))

    mi = minimize_scalar(abs_den, method="bounded", bracket=[0, roots[0]], bounds=[0, roots[0]])
    print(nameU, mi.x)
    if(mi.fun < 1e-3):
        peak_weights[jdx] = np.pi / derivative(denominator, mi.x, dx=1e-6)
        peak_positions[jdx] = mi.x

    #w_vals = 10000
    #w_lin = np.linspace(5e-5, 2, w_vals)
    #ax.plot(w_lin, denominator( w_lin**2 ), color=colors[jdx],
    #            label=f"Termi $@i={np.argmin(deviation_from_inf)+1}$ - $V={nameU}$")
        

ax.plot(names, peak_positions, label="Peak positions")
ax.plot(names, peak_weights, label="Peak weights")

try:
    def func(x, a, b):
        return np.sign(x ) * a * np.abs(x )**b
    
    popt, pcov = curve_fit(func, names, peak_weights, p0=(1, -1))
    v_lin = np.linspace(np.min(names), np.max(names), 200)
    ax.plot(v_lin, func(v_lin, *popt), "--", label=r"$a x^b$")
    ax.set_xlabel(r"$\epsilon$")
    ax.set_ylabel(r"Peak")
    ax.text(0.45, 0.9, f"$a={popt[0]}$", transform = ax.transAxes)
    ax.text(0.45, 0.85, f"$b={popt[1]}$", transform = ax.transAxes)

    popt, pcov = curve_fit(func, names, peak_positions, p0=(1, -1))
    v_lin = np.linspace(np.min(names), np.max(names), 200)
    ax.plot(v_lin, func(v_lin, *popt), "--", label=r"$c x^d$")
    ax.set_xlabel(r"$\epsilon$")
    ax.set_ylabel(r"Peak")
    ax.text(0.45, 0.2, f"$c={popt[0]}$", transform = ax.transAxes)
    ax.text(0.45, 0.15, f"$d={popt[1]}$", transform = ax.transAxes)

except RuntimeError:
    print("Could not estimate curve_fit")

ax.set_xscale("log")
ax.set_yscale("log")

ax.legend()
ax.set_xlabel(r"$V / t$")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
