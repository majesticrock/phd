import matplotlib.pyplot as plt
import numpy as np

M = np.loadtxt(f"data/T0.1/U_modes/-3.00.txt").transpose()
M = 1 / M[abs(M) > 1e-10]

sigma = 0.0005
def gauss(x, mu):
    return np.exp((-(x - mu)**2) / (2*sigma))

off = abs(M[0] - M[int(len(M) / 10)])
lims = [min(M), max(M)] 
lims = [-15, 20]
size = 20000
x_lin = np.linspace(lims[0], lims[1], size)

data = np.zeros(size)
for i in range(0, len(M)):
    data += gauss(x_lin, M[i])

plt.plot(x_lin, data)
#plt.plot(M, "x")

plt.xlim(lims[0], lims[1])
plt.xlabel("$\\epsilon /t$")
plt.show()