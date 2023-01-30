import matplotlib.pyplot as plt
import numpy as np

M = 1/np.loadtxt(f"data/T0.1/U_modes/-3.00.txt").transpose()

sigma = 0.001
def gauss(x, mu):
    return np.exp((-(x - mu)**2) / (2*sigma))

off = abs(M[0] - M[200])
x_lin = np.linspace(np.min(M) - off, np.max(M) + off, 10000)

data = np.zeros(10000)
for i in range(0, len(M)):
    data += gauss(x_lin, M[i])

plt.plot(x_lin, data)
#plt.plot(M, M, "x")

#plt.xlim(-1e5, 1e5)
plt.show()