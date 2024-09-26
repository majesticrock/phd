import numpy as np
import matplotlib.pyplot as plt

def fock(k, k_s):
    return -( 1. + k_s * ( np.atan( (k - 1) / k_s ) - np.atan( (k + 1) / k_s ) ) 
                + (1 + k_s**2 - k**2) / (2 * k) * np.log( ( k_s**2 + (1 + k)**2 ) / ( k_s**2 + (1 - k)**2 ) ) )
    
screenings = np.array([1e-4, 0.1, 0.5, 1., 2.])
k_lin = np.linspace(1e-8, 2., 500)

fig1, ax1 = plt.subplots()

for screening in screenings:
    ax1.plot(k_lin, fock(k_lin, screening), label=f"$k_s / k_\\mathrm{{F}}= {screening}$")

ax1.set_xlabel(r"$k$ $[k_\mathrm{F}]$")
ax1.set_ylabel(r"$\epsilon_\mathrm{Fock}$ $[ \dfrac{e^2 k_\mathrm{F} }{4\pi \epsilon_0} ] $")
ax1.legend()
fig1.tight_layout()


fig2, ax2 = plt.subplots()
screenings = np.linspace(1e-4, 10., 500)
ax2.plot(screenings, fock(1., screenings))
ax2.set_xlabel(r"$k_s / k_\mathrm{F}$")
ax2.set_ylabel(r"$\epsilon_\mathrm{Fock}$ $[ \dfrac{e^2 k_\mathrm{F} }{4\pi \epsilon_0} ] $")
fig2.tight_layout()

plt.show()