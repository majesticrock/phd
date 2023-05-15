import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib as mpl
import sys
import gzip

eps = 1e-5
if(len(sys.argv) > 1):
    data_folder = "data/" + sys.argv[1] + "/"
    name = sys.argv[1]
else:
    name = "T0"
    data_folder = f"data/phases/{name}/"


with gzip.open(data_folder + "cdw_up.dat.gz", 'rt') as f_open:
    CDW_UP = np.loadtxt(f_open)
with gzip.open(data_folder + "cdw_down.dat.gz", 'rt') as f_open:
    CDW_DOWN = np.loadtxt(f_open)
with gzip.open(data_folder + "sc.dat.gz", 'rt') as f_open:
    SC =  abs(np.loadtxt(f_open))
with gzip.open(data_folder + "eta.dat.gz", 'rt') as f_open:
    ETA = abs(np.loadtxt(f_open))

CDW = np.abs(CDW_UP + CDW_DOWN)
AFM = np.abs(CDW_UP - CDW_DOWN)

labels = ["T", "U"]
T_SIZE = len(CDW)
U_SIZE = len(CDW[0])

with gzip.open(data_folder + "cdw_up.dat.gz", 'rt') as fp:
    for i, line in enumerate(fp):
        if i == 2:
            ls = line.split()
            labels[0] = ls[1].split("_")[0]
            U = np.linspace(float(ls[1].split("=")[1]), float(ls[2].split("=")[1]), U_SIZE+1)[:U_SIZE]
        elif i == 3:
            ls = line.split()
            labels[1] = ls[1].split("_")[0]
            T = np.linspace(float(ls[1].split("=")[1]), float(ls[2].split("=")[1]), T_SIZE+1)[:T_SIZE]
        elif i > 3:
            break

for i in range(0, T_SIZE):
    for j in range(0, U_SIZE):
        if(CDW[i][j] > eps):
            CDW[i][j] = 1
        else:
            CDW[i][j] = 0
        if(AFM[i][j] > eps):
            AFM[i][j] = 1
        else:
            AFM[i][j] = 0
        if(SC[i][j] > eps):
            SC[i][j] = 1
        else:
            SC[i][j] = 0
        if(ETA[i][j] > eps):
            ETA[i][j] = 1
        else:
            ETA[i][j] = 0


X, Y = np.meshgrid(U, T)
cmap1 = colors.ListedColormap(['white', 'C0'])
cmap2 = colors.ListedColormap(['white', 'C1'])
cmap3 = colors.ListedColormap(['white', 'C2'])

fig, ax = plt.subplots()

mpl.rcParams["hatch.linewidth"] = 2.5
cset1 = ax.contourf(X, Y, SC, 1, cmap=cmap1, hatches=[None, None])
cset2 = ax.contourf(X, Y, CDW, 1, cmap=cmap2, alpha=0.4)
cset3 = ax.contourf(X, Y, AFM, 1, cmap=cmap3, alpha=0.4)
cset3 = ax.contourf(X, Y, ETA, 1, cmap=cmap2, alpha=0.1)
#cbar = fig.colorbar(cset1)

from matplotlib.patches import Patch

legend_elements = [Patch(facecolor='C0', label=r'SC'),
            Patch(facecolor='C1', label=r'CDW'),
            Patch(facecolor='C2', label=r'AFM')]
ax.legend(handles=legend_elements, loc='upper right')

plt.xlabel(r"$" + labels[0] + "/t$")
plt.ylabel(r"$" + labels[1] + "/t$")

import os
if not os.path.exists("python/build"):
    os.makedirs("python/build")
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}_{name}.pdf")
plt.show()