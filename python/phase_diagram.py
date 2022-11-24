import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib as mpl

eps = 1e-8
data_folder = "data/T/"
CDW = abs(np.loadtxt(data_folder + "cdw.txt"))
SC =  abs(np.loadtxt(data_folder + "sc.txt"))
ETA = abs(np.loadtxt(data_folder + "eta.txt"))

labels = ["T", "U"]
T_SIZE = len(CDW)
U_SIZE = len(CDW[0])

with open(data_folder + "cdw.txt") as fp:
    for i, line in enumerate(fp):
        if i == 2:
            ls = line.split()
            labels[0] = ls[1].split("_")[0]
            U = np.linspace(float(ls[1].split("=")[1]), float(ls[2].split("=")[1]), U_SIZE)
        elif i == 3:
            ls = line.split()
            labels[1] = ls[1].split("_")[0]
            T = np.linspace(float(ls[1].split("=")[1]), float(ls[2].split("=")[1]), T_SIZE)
        elif i > 3:
            break

for i in range(0, T_SIZE):
    for j in range(0, U_SIZE):
        if(CDW[i][j] > eps):
            CDW[i][j] = 1
        else:
            CDW[i][j] = 0
        if(SC[i][j] > eps):
            SC[i][j] = 1
        else:
            SC[i][j] = 0
        if(ETA[i][j] > eps):
            ETA[i][j] = 1
        else:
            ETA[i][j] = 0

X, Y = np.meshgrid(U, T)
cmap1 = colors.ListedColormap(['white', 'green'])
cmap2 = colors.ListedColormap(['white', 'red'])

fig, ax = plt.subplots()

mpl.rcParams["hatch.linewidth"] = 2.5
cset1 = ax.contourf(X, Y-0.5*(T[1]-T[0]), SC, 1, cmap=cmap1, hatches=[None, None])
cset2 = ax.contourf(X, Y-0.5*(T[1]-T[0]), CDW, 1, cmap=cmap2, hatches=[None, r"//"], alpha=0.4)
cset3 = ax.contourf(X, Y-0.5*(T[1]-T[0]), ETA, 1, cmap=cmap2, hatches=[None, r"\\"], alpha=0)
#cbar = fig.colorbar(cset1)
print(U)

from matplotlib.patches import Patch

legend_elements = [Patch(facecolor='green', label=r'$\Delta_{SC}$'),
            Patch(facecolor='red', label=r'$\Delta_{CDW}$')]
ax.legend(handles=legend_elements, loc='upper right')

plt.xlabel(r"$" + labels[0] + "/t$")
plt.ylabel(r"$" + labels[1] + "/t$")
plt.show()