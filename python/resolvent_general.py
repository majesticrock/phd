import numpy as np
import matplotlib.pyplot as plt
import gzip
# Calculates the resolvent in w^2

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

Ts = np.array([0.])
Us = np.array([-2.0])
Vs = np.array([-0.1])

folder = "data/test/"
name_suffix = "phase_sc_"
element_names = ["a", "a+b", "a+ib"]
fig, ax = plt.subplots()

#ax.set_xscale("log")
ax.set_yscale("log")

lss = ["-", "--", "-."]

plot_upper_lim = 8.5

for q, T in enumerate(Ts):
    for r, U in enumerate(Us):
        for s, V in enumerate(Vs):
            name = f"T={T}/U={U}_V={V}/"

            file = f"{folder}{name}one_particle.dat.gz"
            with gzip.open(file, 'rt') as f_open:
                one_particle = np.abs(np.loadtxt(f_open).flatten())

            roots = np.array([np.min(one_particle) * 2, np.max(one_particle) * 2])**2
            a_inf = (roots[0] + roots[1]) * 0.5
            b_inf = ((roots[1] - roots[0]) * 0.25)
            ax.axvspan(np.sqrt(roots[1]), np.sqrt(roots[0]), alpha=.2, color="purple", label="Continuum")

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

            w_vals = 20000
            w_lin = np.linspace(0, plot_upper_lim, w_vals, dtype=complex)**2
            w_lin += 1e-8j
            off = 1

            data = np.zeros(w_vals)

            for idx, element in enumerate(element_names):
                file = f"{folder}{name}resolvent_{name_suffix}{element}.dat.gz"
                with gzip.open(file, 'rt') as f_open:
                    M = np.loadtxt(f_open)
                    A = M[0]
                    B = M[1]

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
                
                if idx == 0:
                    data -= dos( np.copy(w_lin) ).imag
                elif idx == 1:
                    data -= dos( np.copy(w_lin) ).imag
                elif idx == 2:
                    data += dos( np.copy(w_lin) ).imag

            ax.plot(np.sqrt(w_lin.real), data, color=colors[q+r+s],
                linewidth=(plt.rcParams["lines.linewidth"]), label=f"$V={V}$")

legend = plt.legend()

import matplotlib.lines as mlines
dummy_lines = []
dummy_lines.append(mlines.Line2D([],[], color="k", linestyle="-"))
dummy_lines.append(mlines.Line2D([],[], color="k", linestyle="--", linewidth=2*plt.rcParams["lines.linewidth"]))
legend_extra = plt.legend([dummy_lines[i] for i in [0,1]], [r"Amplitude", r"Phase"], loc="upper center")

ax.add_artist(legend)
ax.add_artist(legend_extra)

ax.set_xlabel(r"$\epsilon / t$")
ax.set_ylabel(r"$A(z)$")
fig.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()
