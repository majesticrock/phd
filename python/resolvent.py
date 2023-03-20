import numpy as np
import matplotlib.pyplot as plt

# Calculates the resolvent in 1/w

nameU = "2.00"
folder = "T0"
subfolder = ""

file = f"data/{folder}/V_modes/{subfolder}{nameU}_resolvent.txt"
one_particle = 1 / np.abs(np.loadtxt(f"data/{folder}/V_modes/{subfolder}{nameU}_one_particle.txt").flatten())

M = np.loadtxt(file)
A = M[0]
B = M[1]

w_vals = 10000
w_lin = np.linspace(-15, 15, w_vals, dtype=complex)
w_lin += 1e-1j
off = 1

B_min = 1/16 * ( np.min(one_particle) - np.max(one_particle))**2 #
B_max = 1/16 * ( np.min(one_particle) + np.max(one_particle))**2 #
roots = np.array([np.sqrt((np.sqrt(B_min) - np.sqrt(B_max))**2), np.sqrt((np.sqrt(B_min) + np.sqrt(B_max))**2)])

def r(w):
    ret = np.zeros(len(w), dtype=complex)
    for i in range(0, len(w)):
        root = (w[i]**2 + B_min - B_max)**2 - 4*w[i]**2 * B_min
        if(abs(w[i]) < roots[0]):
            ret[i] = (1/(w[i]*B_min)) * ( w[i]**2 + B_min - B_max + np.sqrt(root, dtype=complex) )
        elif(abs(w[i]) > roots[1]):
            ret[i] = (1/(w[i]*B_min)) * ( w[i]**2 + B_min - B_max - np.sqrt(root, dtype=complex) )
        else:
            ret[i] = (1/(w[i]*B_min)) * ( w[i]**2 + B_min - B_max - np.sqrt(root, dtype=complex) )
    return ret

def dos(w):
    G = w - A[len(A) - off] - B[len(B) - off] #* r( w )
    for j in range(len(A) - off - 1, -1, -1):
        G = w - A[j] - B[j + 1] / G
    return -w * B[0] / G
    
fig, ax = plt.subplots()
ax.plot(w_lin.real, -dos( 1 / w_lin ).imag, "-x", markevery=0.01, label="Imag")
ax.plot(w_lin.real,   -r( 1 / w_lin.real ).imag, label="$r(\\omega)$")
#print(np.trapz(-dos( 1 / w_lin ).imag, dx=20. / w_vals))
#ax.plot(w_lin.real, dos(w_lin).real, "x", label="Real")
#ax.plot(A, 'x', label="$a_i$")
#ax.plot(B, 'o', label="$b_i$")
#ax.set_yscale("log")
#ax.set_ylim(-10, 10)

ax.legend()
ax.set_xlabel(r"$\epsilon / t$")
fig.tight_layout()
plt.show()