import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

def alpha(delta_epsilon):
    return delta_epsilon + 1.

def beta(delta_epsilon):
    return delta_epsilon - 1.

def interaction_CUT(delta_eps, delta_eps_prime):
    A = np.sign(beta(delta_eps_prime))  / ( np.abs(alpha(delta_eps)) + np.abs(beta(delta_eps_prime)) )
    B = np.sign(alpha(delta_eps_prime)) / ( np.abs(beta(delta_eps)) + np.abs(alpha(delta_eps_prime)) )
    return A - B

def interaction_lenz_wegner(delta_eps, delta_eps_prime):
    A = beta(delta_eps_prime)  / ( alpha(delta_eps)**2 + beta(delta_eps_prime)**2  )
    B = alpha(delta_eps_prime) / ( beta(delta_eps)**2  + alpha(delta_eps_prime)**2 )
    return A - B

def interaction_froehlich(delta_eps, irrelevant):
    return 1. / (delta_eps**2 - 1.)

def axis_transform(E, E_P):
    return (E - E_P, E_P)

fig, ax = plt.subplots()
epsilon_space = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(epsilon_space, epsilon_space)

Z = interaction_CUT(*axis_transform(X, Y))
limit = min(-np.min(Z), np.max(Z))
divnorm = colors.TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
contour = ax.contourf(X, Y, Z, levels=1000, cmap='seismic', norm=divnorm)
cbar = fig.colorbar(contour)

ax.set_title("CUT")
ax.set_xlabel(r"$E   (k,  q) / \omega_\mathrm{D}$")
ax.set_ylabel(r"$E_P (k', q) / \omega_\mathrm{D}$")
cbar.set_label(r"$G(k, k', q) / \omega_\mathrm{D}$")

fig2, ax2 = plt.subplots()
Z2 = interaction_lenz_wegner(*axis_transform(X, Y))
limit2 = min(-np.min(Z), np.max(Z))
divnorm2 = colors.TwoSlopeNorm(vmin=-limit2, vcenter=0, vmax=limit2)
contour2 = ax2.contourf(X, Y, Z2, levels=1000, cmap='seismic', norm=divnorm2)
cbar2 = fig2.colorbar(contour2, extend='both')

ax2.set_title("Lenz-Wegner")
ax2.set_xlabel(r"$E   (k,  q) / \omega_\mathrm{D}$")
ax2.set_ylabel(r"$E_P (k', q) / \omega_\mathrm{D}$")
cbar2.set_label(r"$G(k, k', q) / \omega_\mathrm{D}$")

fig3, ax3 = plt.subplots()
Z3 = interaction_froehlich(*axis_transform(X, Y))
limit3 = min(-np.min(Z), np.max(Z))
divnorm3 = colors.TwoSlopeNorm(vmin=-limit3, vcenter=0, vmax=limit3)
contour3 = ax3.contourf(X, Y, Z3, levels=1000, cmap='seismic', norm=divnorm3)
cbar3 = fig3.colorbar(contour3, extend='both')

ax3.set_title("Fröhlich")
ax3.set_xlabel(r"$E   (k,  q) / \omega_\mathrm{D}$")
ax3.set_ylabel(r"$E_P (k', q) / \omega_\mathrm{D}$")
cbar3.set_label(r"$G(k, k', q) / \omega_\mathrm{D}$")

plt.show()