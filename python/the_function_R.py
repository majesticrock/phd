import numpy as np
import matplotlib.pyplot as plt

R = np.loadtxt("data/R.txt")
R_deriv = np.loadtxt("data/R_deriv.txt")
phi = np.linspace(-1, 1, R[0].size, endpoint=True)

fig, axs = plt.subplots(2, 1, sharex=True)

axs[0].plot(phi, R[0], label=r"$\gamma = 0$")
axs[0].plot(phi, R[1], label=r"$\gamma = 0.5$")
axs[0].plot(phi, R[2], label=r"$\gamma = 1$")

axs[1].plot(phi, R_deriv[0], label=r"$\gamma = 0$")
axs[1].plot(phi, R_deriv[1], label=r"$\gamma = 0.5$")
axs[1].plot(phi, R_deriv[2], label=r"$\gamma = 1$")

axs[0].legend()
axs[1].set_xlabel(r"$\varphi$")

axs[0].set_ylabel(r"$R(\varphi, \gamma)$")
axs[1].set_ylabel(r"$\partial_\varphi R(\varphi, \gamma)$")

plt.tight_layout()
plt.show()