import matplotlib.pyplot as plt
import numpy as np

from mrock_centralized_scripts import path_appender
path_appender.append()
import get_data

N=100
df = get_data.load_panda("dwave", "test", "single_gap.json.gz",
                         **get_data.dwave_params(N=N, 
                                                 g=1, 
                                                 V=0.2, 
                                                 E_F=-0.2, 
                                                 omega_D=0.04))

kx = np.linspace(-1, 1, N, endpoint=False)
ky = np.linspace(-1, 1, N, endpoint=False)

X, Y = np.meshgrid(kx, ky)
Z = df["Delta"]
if abs(np.min(Z)) > np.max(Z):
    Z *= -1
Z.shape = (Z.size//N, N)

fig, ax = plt.subplots()
cont = ax.contourf(X, Y, Z, levels=101, cmap='seismic')

cbar = fig.colorbar(cont, ax=ax)
cbar.set_label(r"$\Delta / W$")

ax.set_xlabel(r"$k_x / \pi$")
ax.set_ylabel(r"$k_y / \pi$")

plt.show()