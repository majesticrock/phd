import numpy as np
import matplotlib.pyplot as plt

V = np.array([-0.05, -0.1, -0.15, -0.2, -0.25, -0.3, -0.35])
d300 = np.array([-1.257060e-03, -8.089023e-03, -1.497539e-02, -2.192384e-02, -2.894403e-02, -3.604763e-02, -4.324850e-02])
d900 = np.array([-3.074002e-04, -2.688230e-04, -2.326862e-04, -1.869452e-03, -6.207933e-03, -1.088450e-02, -1.592057e-02])

plt.plot(V, d300, label="300")
plt.plot(V, d900, label="900")
plt.legend()
plt.show()