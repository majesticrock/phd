import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib as mpl

eps = 1e-6

CDW = (np.loadtxt("data/basic_hubbard_cdw.txt"))
SC =  (np.loadtxt("data/basic_hubbard_sc.txt"))
ETA = (np.loadtxt("data/basic_hubbard_eta.txt"))

T = np.loadtxt("data/basic_hubbard_T.txt")
U = np.loadtxt("data/basic_hubbard_U.txt")

T_SIZE = len(CDW)
U_SIZE = len(CDW[0])
Z = np.zeros((T_SIZE, U_SIZE))

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
cset1 = ax.contourf(X-0.5*(U[1]-U[0]), Y-0.5*(T[1]-T[0]), SC, 1, cmap=cmap1, hatches=[None, None])
cset2 = ax.contourf(X-0.5*(U[1]-U[0]), Y-0.5*(T[1]-T[0]), CDW, 1, cmap=cmap2, hatches=[None, r"//"], alpha=0.4)
cset3 = ax.contourf(X-0.5*(U[1]-U[0]), Y-0.5*(T[1]-T[0]), ETA, 1, cmap=cmap2, hatches=[None, r"\\"], alpha=0)
#cbar = fig.colorbar(cset1)

from matplotlib.patches import Patch

legend_elements = [Patch(facecolor='green', label=r'$\Delta_{SC}$'),
            Patch(facecolor='red', label=r'$\Delta_{CDW}$')]
ax.legend(handles=legend_elements, loc='upper right')

plt.xlabel(r"$U / t$")
plt.ylabel(r"$V / t$")
plt.show()