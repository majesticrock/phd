import matplotlib.pyplot as plt
import numpy as np

M = np.loadtxt(f"data/T0.1/U_modes/-100.00.txt").transpose()
E = np.loadtxt(f"data/T0.1/U_modes/-100.00_one_particle.txt").transpose()
k_lin = np.linspace(-1, 1, len(M[0]), endpoint=False)
for i in range(0, len(M)):
    plt.plot(k_lin, M[i], label=i)
#for i in range(0, len(E)):
#    plt.plot(k_lin, 2*E[i], "--")

#x = np.append(k_lin, 2 + k_lin)
#y = np.append(M[6], M[6])
#plt.plot(x, y)
#from scipy.optimize import curve_fit
#def f(x, a, b, c, d):
#    return np.cos(np.pi*a*x + d) * b + c
#
#popt, pcov = curve_fit(f, x, y)
#plt.plot(x, f(x, *popt), "--")
#print(popt)
#
def e0(x, y):
    return 2*(np.cos(np.pi*x) + np.cos(np.pi*y))
#plt.plot(k_lin, 4*e0(k_lin, 0), "k:")

plt.xlabel(r"$k/\pi$")
plt.ylabel(r"$\epsilon / t$")
plt.tight_layout()

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.svg")
plt.show()
