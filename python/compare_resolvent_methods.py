import numpy as np
import matplotlib.pyplot as plt

# Calculates the resolvent in 1/w

nameU = "-0.10"
folder = "T0"
subfolders = ["_L", "_BL"]
labels = ["$L=30$", "$L=50$"]
lss = ["-", "--", ":"]

fig, ax = plt.subplots()
for s, subfolder in enumerate(subfolders):
    file = f"data/{folder}/V_modes/{nameU}_resolvent{subfolder}.txt"
    one_particle = 1 / np.abs(np.loadtxt(f"data/{folder}/V_modes/{nameU}_one_particle{subfolder}.txt").flatten())

    M = np.loadtxt(file)
    A = M[0]
    B = M[1]

    w_vals = 10000
    w_lin = 1 / np.linspace(-10, 10, w_vals, dtype=complex)
    w_lin += 5e-2j
    off = 1

    def dos(w):
        G = w - A[len(A) - off] - B[len(B) - off]
        for j in range(len(A) - off - 1, -1, -1):
            G = w - A[j] - B[j + 1] / G
        return w * B[0] / G

    ax.plot(1 / w_lin.real, -dos( w_lin ).imag, "-", label=f"{labels[s]} Lanczos")
    R = np.loadtxt(f"data/{folder}/V_modes/{nameU}{subfolder}.txt")
    ax.plot(np.linspace(-10, 10, len(R)), R, "--", label=f"{labels[s]} Exact")

B_min = 1/16 * ( np.min(one_particle) - np.max(one_particle))**2 #
B_max = 1/16 * ( np.min(one_particle) + np.max(one_particle))**2 #
roots = np.array([np.sqrt((np.sqrt(B_min) - np.sqrt(B_max))**2), np.sqrt((np.sqrt(B_min) + np.sqrt(B_max))**2)])
ax.axvspan(-1/roots[1], -1/roots[0], alpha=.2, color="purple", label="Continuum")
ax.axvspan(1/roots[1], 1/roots[0], alpha=.2, color="purple")

ax.legend()
ax.set_xlabel(r"$\epsilon / t$")
fig.tight_layout()
plt.show()