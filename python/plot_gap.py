import numpy as np
import matplotlib.pyplot as plt

CDW = (np.loadtxt("data/basic_hubbard_cdw.txt"))
SC =  (np.loadtxt("data/basic_hubbard_sc.txt"))
ETA = (np.loadtxt("data/basic_hubbard_eta.txt"))

T_SIZE = len(CDW)
U_SIZE = len(CDW[0])

with open("data/basic_hubbard_cdw.txt") as fp:
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

fig, ax = plt.subplots(2, 1)

for i in range(0, 20, 2):
    ax[0].plot(T, data.transpose()[i], label="$U=" + str(U[i]) + "$")

for i in range(0, 20, 2):
    ax[1].plot(U, data[i], label="$T=" + str(T[i]) + "$")


ax[0].set_xlabel(r"$T/t$")
ax[1].set_xlabel(r"$U/t$")
ax[0].set_ylabel(r"$\Delta_{tot}$")
ax[1].set_ylabel(r"$\Delta_{tot}$")
#ax[0].legend()
#ax[1].legend()
plt.tight_layout()
plt.show()