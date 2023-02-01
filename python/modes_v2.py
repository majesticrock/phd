import matplotlib.pyplot as plt
import numpy as np

M = np.loadtxt(f"data/T0.1/U_modes/0.00.txt").transpose()
M = 1 / M[abs(M) > 1e-6]

sigma = 0.01
def gauss(x, mu):
    return np.exp((-(x - mu)**2) / (2*sigma))

lims = [min(M) - 1, max(M) + 1]
#lims = [-450, 450]
size = 20000
x_lin = np.linspace(lims[0], lims[1], size)

data = np.zeros(size)
for i in range(0, len(M)):
    data += gauss(x_lin, M[i])

plt.plot(x_lin, data)
#plt.plot(M, "x")

#plt.xlim(lims[0], lims[1])
plt.xlabel("$\\epsilon /t$")
plt.show()