import numpy as np
import matplotlib.pyplot as plt

# Calculates the resolvent in 1/w

file = "data/resolvent.txt"

M = np.loadtxt(file)
A = M[0]
B = M[1]

w_vals = 2000
w_lin = np.linspace(-12, 12, w_vals, dtype=complex)
w_lin += 1e-3j
off = 1

B_min = B[len(B) - off]
B_max = B[len(B) - off - 1]

roots = [np.sqrt((np.sqrt(B_min) - np.sqrt(B_max))**2), np.sqrt((np.sqrt(B_min) + np.sqrt(B_max))**2)]

def r(w):
    #w = wp.real
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

def denominator(w):
    G = 1. - w * A[len(A) - off] - w * B[len(B) - off] * r( 1/w )
    for j in range(len(A) - off - 1, -1, -1):
        G = 1. - w * A[j] - w**2 * B[j + 1] / G
    return G

def dos(w):
    G = 1. - w * A[len(A) - off] - w * B[len(B) - off] * r( 1/w )
    for j in range(len(A) - off - 1, -1, -1):
        G = 1. - w * A[j] - w**2 * B[j + 1] / G
    return -B[0] / G
    
fig, ax = plt.subplots()
ax.plot(w_lin.real, -dos(w_lin).imag)
ax.plot(w_lin.real, -r(1/w_lin).imag)
#ax.plot(w_lin.real, dos(w_lin).real, "x")
#ax.plot(A[:100], linestyle="-", marker='x')
#ax.plot(B[:100], linestyle="-", marker='o')


ax.set_xlabel(r"$\epsilon / t$")
fig.tight_layout()
plt.show()