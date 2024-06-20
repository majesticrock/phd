import numpy as np
import matplotlib.pyplot as plt

# Define the function
def f(x, a, n):
    return np.sinh(2 * (a / n) * x - a ) / np.sinh(a) + 1

# Define n
n = 1000

# Define the x range
x = np.linspace(-0.1 * n, 1.1 * n, 400)

# Plotting for different values of a
plt.figure(figsize=(10, 8))

# Plotting for different values of a
for a in [0.5, 1.0, 2.0, 4.0, 10., 100.]:
    y = f(x, a, n)
    plt.plot(x, y, label=f'a = {a:.1f}')

# Plotting the fixed points
plt.plot([0, 0.5 * n, n], [0, 1, 2], 'ro')  # Points (-n, -1), (0, 0), (n, 1)

plt.title('Function f(x) with Adjustable Steepness')
plt.xlabel('x')
plt.ylabel('y')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.legend()
plt.grid(True)
plt.ylim(-0.5, 2.5)

plt.show()
