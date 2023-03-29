import numpy as np
import matplotlib.pyplot as plt

# Calculates the resolvent in 1/w

nameU = "0.10"
folder = "T0_test"
subfolders = ["ev", "ldlt", "llt"]
lss = ["-", "--", ":"]

fig, ax = plt.subplots()
for s, subfolder in enumerate(subfolders):
    file = f"data/{folder}/{subfolder}/{nameU}_resolvent.txt"
    one_particle = 1 / np.abs(np.loadtxt(f"data/{folder}/{subfolder}/{nameU}_one_particle.txt").flatten())

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

    ax.plot(1 / w_lin.real, -dos( w_lin ).imag, "-", label=f"{subfolder} Lanczos", linestyle=lss[s])
    R = np.loadtxt(f"data/{folder}/{subfolder}/{nameU}.txt")
    ax.plot(np.linspace(-10, 10, len(R)), R, "--", label=f"{subfolder} Exact", linestyle=lss[s])

ax.legend()
ax.set_xlabel(r"$\epsilon / t$")
fig.tight_layout()
plt.show()