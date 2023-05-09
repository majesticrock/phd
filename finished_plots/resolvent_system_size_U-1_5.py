import numpy as np
import matplotlib.pyplot as plt
import gzip
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([-1.5])
Vs = np.array([-0.5])

folders = ["L=30", "L=40", "L=50", "L=60", "L=70"]
name_suffix = "SC"
fig, ax = plt.subplots()
types = [ "phase"]#"higgs",
lss = ["-", "-", "-", "--", ":", ]
plt.rcParams["lines.linewidth"] *= 2

ax.set_yscale("log")
axins = ax.inset_axes([0.1, 0.04, 0.2, 0.4])

for folder, ls in zip(folders, lss):
    for q, T in enumerate(Ts):
        for r, U in enumerate(Us):
            for s, V in enumerate(Vs):
                name = f"T={T}/U={U}_V={V}/"

                file = f"data/{folder}/{name}one_particle.dat.gz"
                with gzip.open(file, 'rt') as f_open:
                    one_particle = np.abs(np.loadtxt(f_open).flatten())

                print("Gap = ", np.min(one_particle))
                roots = np.array([np.min(one_particle) * 2, np.max(one_particle) * 2])**2
                a_inf = (roots[0] + roots[1]) * 0.5
                b_inf = ((roots[1] - roots[0]) * 0.25)

                if(ls==":"): ax.axvspan(np.sqrt(roots[1]), np.sqrt(roots[0]), alpha=.2, color="purple", label="Continuum")
                for idx, type in enumerate(types):
                    file = f"data/{folder}/{name}resolvent_{type}_{name_suffix}.dat.gz"
                    with gzip.open(file, 'rt') as f_open:
                        M = np.loadtxt(f_open)
                        A = M[0][:20]
                        B = M[1][:20]

                    w_vals = 20000
                    w_lin = np.linspace(0, 9, w_vals, dtype=complex)**2
                    w_lin += 1e-3j
                    off = 1

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

                    off_termi = len(A) - 1 - np.argmin(deviation_from_inf)
                    #print("Terminating at i=", np.argmin(deviation_from_inf))
                    def dos(w):
                        for i in range(0, len(w)):
                            if(w[i].real > roots[0] and w[i].real < roots[1]):
                                w[i] = w[i].real

                        G = w - A[len(A) - off_termi] - B[len(B) - off_termi] * terminator( w )
                        for j in range(len(A) - off_termi - 1, -1, -1):
                            G = w - A[j] - B[j + 1] / G
                        return B[0] / G

                    ax.plot(np.sqrt(w_lin.real), -dos( np.copy(w_lin) ).imag, 
                            linestyle=ls, label=f"${folder}$")
                    axins.plot(np.sqrt(w_lin.real), -dos( np.copy(w_lin) ).imag, 
                            linestyle=ls, label=f"${folder}$")


ax.legend()

ax.set_xlabel(r"$\epsilon / t$")
ax.set_ylabel(r"$A(z)$")

x1, x2, y1, y2 = 0, 0.25, 0.63, 440
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticklabels([])
axins.set_yticklabels([])
ax.indicate_inset_zoom(axins, edgecolor="black")

fig.tight_layout()

import os
plt.savefig(f"finished_plots/build/{os.path.basename(__file__).split('.')[0]}.svg")
