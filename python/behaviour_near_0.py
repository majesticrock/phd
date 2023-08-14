import numpy as np
import matplotlib.pyplot as plt
import sys
import gzip

if(len(sys.argv) > 1):
    data_folder = "data/" + sys.argv[1] + "/"
    name = sys.argv[1]
else:
    name = "T0"
    data_folder = f"data/phases/square/{name}/"

with gzip.open(data_folder + "cdw.dat.gz", 'rt') as f_open:
    CDW = abs(np.loadtxt(f_open))
with gzip.open(data_folder + "sc.dat.gz", 'rt') as f_open:
    SC =  abs(np.loadtxt(f_open))
with gzip.open(data_folder + "afm.dat.gz", 'rt') as f_open:
    AFM = abs(np.loadtxt(f_open))

labels = ["T", "U"]
T_SIZE = len(CDW)

with gzip.open(data_folder + "cdw.dat.gz", 'rt') as fp:
    for i, line in enumerate(fp):
        if i == 3:
            ls = line.split()
            labels[0] = ls[1].split("_")[0]
            T = np.linspace(float(ls[1].split("=")[1]), float(ls[2].split("=")[1]), T_SIZE+1)[:T_SIZE]
        elif i > 3:
            break

data = np.sqrt(CDW*CDW + SC*SC + AFM*AFM)
data = data.transpose()

plt.plot(T, np.log(data), label='Mean Field')

def theory(u, a):
    u = np.abs(u)
    return np.log(a * (4. / u) * np.exp(-2 * np.pi * np.sqrt(1. / u)))

from scipy.optimize import curve_fit
#popt, pcov = curve_fit(theory, T, data)
plt.plot(T, theory(T, 1), label="Kopietz")

plt.xlabel('$' + labels[0] + '/t$')
plt.ylabel(r'$\Delta/t$')
plt.legend()
plt.tight_layout()
plt.show()