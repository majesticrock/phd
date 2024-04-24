import numpy as np
import matplotlib.pyplot as plt

R = np.loadtxt("data/integration_by_parts.txt")
gamma = np.linspace(-3, 3, R[0].size, endpoint=True)

plt.plot(gamma, R[0], label=r"$I_1$")
plt.plot(gamma, -0.5 * R[1], label=r"$-I_2$")
plt.plot(gamma, R[0] - 0.5 * R[1], label=r"$I_1 - I_2$")


plt.legend()
plt.xlabel(r"$\gamma$")
plt.ylabel(r"$I(\gamma)$")

plt.tight_layout()
plt.show()