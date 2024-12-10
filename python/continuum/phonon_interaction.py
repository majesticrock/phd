import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

LEVELS = 400
N_EPSILON = 100

X_LABEL = r"$\Delta \epsilon (k, q) / \omega_\mathrm{D}$"
Y_LABEL = r"$\Delta \epsilon (k', q) / \omega_\mathrm{D}$"
G_LABEL = r"$G(k, k', q) / (M^2 / \omega_\mathrm{D})$"

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

def symmetrize_interaction(interaction, delta_eps, delta_eps_prime):
    return 0.5 * (interaction(delta_eps, delta_eps_prime) + interaction(delta_eps_prime, delta_eps))
    
def axis_transform(E, E_P):
    return (E, E_P)

fig, ax = plt.subplots()
epsilon_space = np.linspace(-5, 5, N_EPSILON)
X, Y = np.meshgrid(epsilon_space, epsilon_space)

Z = symmetrize_interaction(interaction_CUT, *axis_transform(X, Y))
limit = min(-np.min(Z), np.max(Z))
divnorm = colors.TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
contour = ax.contourf(X, Y, Z, levels=LEVELS, cmap='seismic', norm=divnorm)
cbar = fig.colorbar(contour)

ax.set_title("CUT")
ax.set_xlabel(X_LABEL)
ax.set_ylabel(Y_LABEL)
cbar.set_label(G_LABEL)

fig2, ax2 = plt.subplots()
Z2 = symmetrize_interaction(interaction_lenz_wegner, *axis_transform(X, Y))
limit2 = min(-np.min(Z), np.max(Z))
divnorm2 = colors.TwoSlopeNorm(vmin=-limit2, vcenter=0, vmax=limit2)
contour2 = ax2.contourf(X, Y, Z2, levels=LEVELS, cmap='seismic', norm=divnorm2)
cbar2 = fig2.colorbar(contour2, extend='both')

ax2.set_title("Lenz-Wegner")
ax2.set_xlabel(X_LABEL)
ax2.set_ylabel(Y_LABEL)
cbar2.set_label(G_LABEL)

fig3, ax3 = plt.subplots()
Z3 = symmetrize_interaction(interaction_froehlich, *axis_transform(X, Y))
limit3 = min(-np.min(Z), np.max(Z))
divnorm3 = colors.TwoSlopeNorm(vmin=-limit3, vcenter=0, vmax=limit3)
contour3 = ax3.contourf(X, Y, Z3, levels=LEVELS, cmap='seismic', norm=divnorm3)
cbar3 = fig3.colorbar(contour3, extend='both')

ax3.set_title("Fröhlich")
ax3.set_xlabel(X_LABEL)
ax3.set_ylabel(Y_LABEL)
cbar3.set_label(G_LABEL)

plt.show()