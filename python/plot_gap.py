import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap

data_folder = "data/V/"

CDW = (np.loadtxt(data_folder + "cdw.txt"))
SC =  (np.loadtxt(data_folder + "sc.txt"))
ETA = (np.loadtxt(data_folder + "eta.txt"))

T_SIZE = len(CDW)
U_SIZE = len(CDW[0])

with open(data_folder + "cdw.txt") as fp:
    for i, line in enumerate(fp):
        if i == 2:
            ls = line.split()
            U = np.linspace(float(ls[1].split("=")[1]), float(ls[2].split("=")[1]), U_SIZE)
        elif i == 3:
            ls = line.split()
            T = np.linspace(float(ls[1].split("=")[1]), float(ls[2].split("=")[1]), T_SIZE)
        elif i > 3:
            break

data = np.sqrt(CDW*CDW + SC*SC + ETA*ETA)

length = 32
colors = cm.gist_rainbow(np.linspace(0, 1, length))
my_cmap = ListedColormap(colors[:,:-1])
fig, ax = plt.subplots(1, 1, figsize=(9,6), constrained_layout=True)
sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=U[0], vmax=U[length - 1]))


for i in range(0, length, 1):
    ax.plot(T, data.transpose()[i], c=colors[i], label="$U=" + str(U[i]) + "$")

cbar = fig.colorbar(sm, ax=ax, label='$U/t$')
cbar.set_ticks(U[0::8])
cbar.set_ticklabels(U[0::8])

ax.set_xlabel(r"$T/t$")
ax.set_ylabel(r"$\Delta_{tot}$")

plt.show()