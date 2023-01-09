from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


k = np.linspace(-np.pi, np.pi, 200)
Z = np.zeros((len(k), len(k)))
for i in range(0, len(k)):
    for j in range(0, len(k)):
        Z[i][j] = -2*np.cos(k[i]) -2*np.cos(k[j])

X, Y = np.meshgrid(k, k)

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

plt.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()