import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import colors
import sys
import gzip

eps = 1e-10
if(len(sys.argv) > 1):
    data_folder = "data/" + sys.argv[1] + "/"
    name = sys.argv[1]
else:
    name = "V+2"#"T0_L200"
    data_folder = f"data/phases/square/{name}/"

swapAxis = False

def pair_sort(pair_arr, sortBy):
    if len(pair_arr) < 2: 
        return
    n = len(pair_arr[sortBy])
    for i in range(sortBy, n):
        for j in range(i + 1, n):
            if pair_arr[sortBy][i] > pair_arr[sortBy][j]:
                pair_arr[sortBy][i], pair_arr[sortBy][j] = pair_arr[sortBy][j], pair_arr[sortBy][i]
                pair_arr[1 - sortBy][i], pair_arr[1 - sortBy][j] = pair_arr[1 - sortBy][j], pair_arr[1 - sortBy][i]

file_names = np.array(["cdw", "afm", "sc", "xi_sc"])
crudeData = []

for fname in file_names:
    with gzip.open(data_folder + f"{fname}.dat.gz", 'rt') as f_open:
        if swapAxis:
            crudeData.append(np.loadtxt(f_open).transpose())
        else:
            crudeData.append(np.loadtxt(f_open))

labels = ["T", "U"]
T_SIZE = len(crudeData[0])
U_SIZE = len(crudeData[0][0])

with gzip.open(data_folder + "cdw.dat.gz", 'rt') as fp:
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

for k in range(0, len(file_names)):
    for i in range(0, T_SIZE):
        for j in range(0, U_SIZE):
            if(crudeData[k][i][j] > eps):
                crudeData[k][i][j] = 1
            else:
                crudeData[k][i][j] = 0


X, Y = np.meshgrid(U, T)
cmaps = []
cmaps.append(colors.ListedColormap([colors.to_rgba('white', 0), colors.to_rgba('C0', 0.5)]))
cmaps.append(colors.ListedColormap([colors.to_rgba('white', 0), colors.to_rgba('C1', 0.5)]))
cmaps.append(colors.ListedColormap([colors.to_rgba('white', 0), colors.to_rgba('C2', 0.5)]))
cmaps.append(colors.ListedColormap([colors.to_rgba('white', 0), colors.to_rgba('C3', 0.5)]))

fig, ax = plt.subplots()

for i in range(0, len(file_names)):
    if swapAxis:
        ax.contourf(X, Y, crudeData[i], 1, cmap=cmaps[i])
    else:
        ax.contourf(Y, X, crudeData[i], 1, cmap=cmaps[i])

from matplotlib.patches import Patch

legend_elements = [Patch(facecolor='C0', label=r'CDW'),
            Patch(facecolor='C1', label=r'AFM'),
            Patch(facecolor='C2', label=r'$s$')]
            #Patch(facecolor='C3', label=r'$d_{x^2 - y^2}$'),
            #Line2D([0], [0], label='Micnas', color='k', linestyle="--")]
            #,Patch(facecolor='C4', label=r'$\tilde{s}$')]
ax.legend(handles=legend_elements, loc='upper left')

micnas_d   =  np.loadtxt("data/micnas_d_wave.csv").transpose()
micnas_afm = np.loadtxt("data/micnas_cdw_afm.csv").transpose()

#ax.plot(micnas_d[0], micnas_d[1], "k--")
#ax.plot(micnas_afm[0], micnas_afm[1], "k--", label="Micnas")
#
#ax.set_xlim(-2, 2)
#ax.set_ylim(-2, 2)

plt.xlabel(r"$" + labels[1] + "/t$")
plt.ylabel(r"$" + labels[0] + "/t$")

import os
if not os.path.exists("python/build"):
    os.makedirs("python/build")
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}_{name}.pdf")
plt.show()