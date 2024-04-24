import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
import gzip

EPS = 1e-8

data_folder = "data/modes/square/dos_64k/T=0.0/U=-2.0/V=-0.1"
filename = "resolvent_phase_SC"

file = f"{data_folder}/one_particle.dat.gz"
with gzip.open(file, 'rt') as f_open:
    one_particle = np.abs(np.loadtxt(f_open).flatten())
    
roots = np.array([np.min(one_particle) * 2, np.max(one_particle) * 2])**2
a_infinity = (roots[0] + roots[1]) * 0.5
b_infinity = (roots[1] - roots[0]) * 0.25

file = f"{data_folder}/{filename}.dat.gz"
with gzip.open(file, 'rt') as f_open:
    M = np.loadtxt(f_open)
    A = M[0]
    B = M[1]
    
deviation_from_infinity = np.zeros(len(A) - 1)
for i in range(0, len(A) - 1):
    deviation_from_infinity[i] = abs((A[i] - a_infinity) / a_infinity) + abs((np.sqrt(B[i + 1]) - b_infinity) / b_infinity)
terminate_at = len(A) - np.argmin(deviation_from_infinity)
print("Terminating at i =", np.argmin(deviation_from_infinity))

A = unp.uarray(M[0], np.full(len(M[0]), EPS))
B = unp.uarray(M[1], np.full(len(M[1]), EPS))

def continued_fraction(w_param):
    w = w_param**2
    G = w - A[len(A) - terminate_at]
    for j in range(len(A) - terminate_at - 1, -1, -1):
        G = w - A[j] - B[j + 1] / G
    return B[0] / G

def diverges_when(w_param):
    w = w_param**2
    G = w - A[len(A) - terminate_at]
    for j in range(len(A) - terminate_at - 1, 0, -1):
        G = w - A[j] - B[j + 1] / G
    return A[0] - B[1] / G

num_values = 5000
w_space = np.linspace(0, 9, num_values)
uw_space = unp.uarray(w_space, np.full(num_values, 0))

res = diverges_when(uw_space)
std_dev = np.zeros(num_values)
values = np.zeros(num_values)
for i in range(0, num_values):
    values[i] = res[i].nominal_value
    std_dev[i] = res[i].std_dev if np.abs(values[i] - w_space[i]**2) < 1 else 0

import matplotlib.pyplot as plt
plt.plot(w_space, std_dev)

plt.axvspan(np.sqrt(roots[0]), np.sqrt(roots[1]), alpha=.2, color="purple", label="Continuum")
#plt.yscale("log")
plt.xlabel("$z / t$")
plt.ylabel("Uncertainty")
plt.show()