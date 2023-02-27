from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


k = np.linspace(-1, 1, 200)
Z = np.zeros((len(k), len(k)))
for i in range(0, len(k)):
    for j in range(0, len(k)):
        Z[i][j] = -2*np.cos(np.pi*k[i]) -2*np.cos(np.pi*k[j])

X, Y = np.meshgrid(k, k)

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, antialiased=False, edgecolor='none', linewidth=0, zorder=1)
ax.plot(-k[100:], k[:100], np.zeros(100), color="black", zorder=20000, linewidth=2)
ax.plot(k[100:],  k[:100], np.zeros(100), color="black", zorder=20000, linewidth=2)
ax.plot(k[:100],  k[100:], np.zeros(100), color="black", zorder=20000, linewidth=2)
ax.plot(-k[:100], k[100:], np.zeros(100), color="black", zorder=20000, linewidth=2)

ax.set_xlabel("$k_x / \\pi$")
ax.set_ylabel("$k_y / \\pi$")
ax.set_zlabel("$\epsilon / t$")

plt.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()