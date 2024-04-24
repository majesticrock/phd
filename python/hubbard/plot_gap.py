import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
import sys
import gzip

if(len(sys.argv) > 1):
    data_folder = "data/" + sys.argv[1] + "/"
    name = sys.argv[1]
else:
    name = "dos_T0"
    data_folder = f"data/phases/square/{name}/"

x_axis_is_first = True

with gzip.open(data_folder + "cdw.dat.gz", 'rt') as f_open:
    CDW = abs(np.loadtxt(f_open))
with gzip.open(data_folder + "sc.dat.gz", 'rt') as f_open:
    SC =  abs(np.loadtxt(f_open))
with gzip.open(data_folder + "afm.dat.gz", 'rt') as f_open:
    AFM = abs(np.loadtxt(f_open))

labels = ["T", "U"]
T_SIZE = len(CDW)
U_SIZE = len(CDW[0])
print(T_SIZE, U_SIZE)

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

data = np.sqrt(CDW*CDW + SC*SC + AFM*AFM)

if(not x_axis_is_first):
    T_SIZE, U_SIZE = U_SIZE, T_SIZE
    U, T = T, U
    data = data.transpose()
    labels[0], labels[1] = labels[1], labels[0]

length = len(T)
colors = cm.gist_rainbow(np.linspace(0, 1, length))
my_cmap = ListedColormap(colors[:,:-1])
fig, ax = plt.subplots(1, 1, constrained_layout=True)
sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=T[0], vmax=T[length - 1]))


for i in range(0, length, 1):
    ax.plot(U, data[i], c=colors[i])

cbar = fig.colorbar(sm, ax=ax, label='$' + labels[1] + '/t$')
cbar_ticks = np.round(np.linspace(T[0], T[-1] - (T[0] - T[1]), 11), 4)
cbar.set_ticks(cbar_ticks)
cbar.set_ticklabels(cbar_ticks)

ax.set_xlabel(r"$" + labels[0] + "/t$")
ax.set_ylabel(r"$\Delta_{tot}$")

import os
if not os.path.exists("python/build"):
    os.makedirs("python/build")
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}_{name}.svg")
plt.show()