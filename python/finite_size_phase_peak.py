import numpy as np
import matplotlib.pyplot as plt

N = np.array([300, 600, 900, 1200, 1500, 3000])
W = np.array([0.03163, 0.02114, 0.017246, 0.015029, 0.01357, 0.00969])

plt.plot(1/N, W, "o-")
plt.show()