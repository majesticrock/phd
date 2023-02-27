import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

nameU = "-2.00"
N = np.loadtxt(f"data/T0.1/U_modes/{nameU}_one_particle.txt")
M_low = N[:len(N[0]/2)]
M_high = N[len(N[0]/2):]

k = np.linspace(-np.pi, np.pi, len(M_low[0]), endpoint=False)
X, Y = np.meshgrid(k, k)
fig = plt.figure()
ax = fig.gca(projection='3d')
surf  = ax.plot_surface(X, Y, M_low,  cmap=cm.viridis, vmin=np.min(M_low), vmax=np.max(M_high), antialiased=False, edgecolor='none', linewidth=0)
surf2 = ax.plot_surface(X, Y, M_high, cmap=cm.viridis, vmin=np.min(M_low), vmax=np.max(M_high), antialiased=False, edgecolor='none', linewidth=0)

ax.set_xlabel("$k_x / \\pi$")
ax.set_ylabel("$k_y / \\pi$")
ax.set_zlabel("$\epsilon / t$")

plt.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()