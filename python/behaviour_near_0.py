import numpy as np
import matplotlib.pyplot as plt
import gzip

data_folder = "data/phases/square/T0/afm.dat.gz"#f"data/phases/small_U/afm_square.dat.gz"

with gzip.open(data_folder, 'rt') as f_open:
    AFM = abs(np.loadtxt(f_open))

labels = ["T", "U"]
T_SIZE = len(AFM)

with gzip.open(data_folder, 'rt') as fp:
    for i, line in enumerate(fp):
        if i == 3:
            ls = line.split()
            labels[0] = ls[1].split("_")[0]
            T = np.linspace(float(ls[1].split("=")[1]), float(ls[2].split("=")[1]), T_SIZE+1)[:T_SIZE]
        elif i > 3:
            break

AFM = AFM.transpose()

plt.plot(T, np.log(AFM), label='Mean Field - Square')

def theory(u, a):
    u = np.abs(u)
    return np.log(a * (4. / u) * np.exp(-2 * np.pi * np.sqrt(1. / u)))
plt.plot(T, theory(T, 1), "--", label="Kopietz")

data_folder = f"data/phases/small_U/afm_cube.dat.gz"

with gzip.open(data_folder, 'rt') as f_open:
    AFM = abs(np.loadtxt(f_open))

labels = ["T", "U"]
AFM = AFM.transpose()

plt.plot(T, np.log(AFM), label='Mean Field - SC')

plt.xlabel('$' + labels[1] + '/t$')
plt.ylabel(r'$\ln(\Delta/t)$')
plt.legend()
plt.tight_layout()

import os
if not os.path.exists("python/build"):
    os.makedirs("python/build")
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()