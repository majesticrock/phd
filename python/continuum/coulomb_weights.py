import numpy as np
import matplotlib.pyplot as plt

def f(q, k, lam):
    return q * np.log((lam**2 + (q + k)**2) / (lam**2 + (q - k)**2)) / k

x = np.linspace(0., 2., 5000)
plt.figure()

for k in [0.01, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75]:
    y = f(x, k, 1e-3)
    plt.plot(x, y, label=f'$k = {k:.2f} k_\\mathrm{{F}}$')


plt.xlabel('$q / k_\\mathrm{{F}}$')
plt.ylabel('$w(k,q)$')
#plt.legend(ncol=3, loc="upper center")
plt.grid(True)

plt.show()
