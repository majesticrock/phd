import numpy as np
import matplotlib.pyplot as plt
import os


import mrock.continued_fraction as cf
from mrock.get_data import *
# Create a DataLoader instance.
# By default, this will look for data in the package's default data directory
# or in the directory specified by the MROCK_DATA_DIR environment variable.
data_loader = DataLoader()

# Load the superconducting test case.
#
# This uses the cube geometry/data set and a small attractive
# non-local interaction V.
pd_data = data_loader.load_panda(
    model="hubbard",
    subdir=os.path.join("cube", "test"),
    file="resolvents.json.gz",
    **hubbard_params(0.0, -2.5, -0.1)
)


# Create a ContinuedFraction evaluator.
#
# The object reads the continued-fraction coefficients from pd_data and
# uses them to evaluate the corresponding resolvents.
resolvents = cf.ContinuedFraction(pd_data)

fig, ax = plt.subplots()
ax.set_ylim(-0.05, 1.)
ax.set_xlabel(r"$\omega / t$")
ax.set_ylabel(r"$\mathcal{A} (\omega)  / t^{-1}$")
colors = np.array((
    "blue",
    "orange",
    "black",
    "limegreen",
    "deepskyblue",
    "magenta"
))
linestyles = ["-", "-.", "--", "-", "--", ":"]

# Construct the frequency grid at which the spectral functions are
# evaluated.
w_lin = np.linspace(
    -0.01,
    pd_data["continuum_boundaries"][1] + 0.3,
    5000,
    dtype=complex
)
w_lin += 1e-6j

# Plot the superconducting phase spectral density.
#
# The resolvent name "phase_SC" refers to the superconducting phase channel.
ax.plot(
    w_lin.real,
    resolvents.spectral_density(w_lin, "phase_SC"),
    label="Phase",
    ls=linestyles[0],
    c=colors[0]
)

# Plot the superconducting amplitude, or Higgs, spectral density.
# The resolvent name "amplitude_SC" refers to the amplitude channel of the superconducting order parameter.
ax.plot(
    w_lin.real,
    resolvents.spectral_density(w_lin, "amplitude_SC"),
    label="Higgs",
    ls=linestyles[1],
    c=colors[1]
)

# Plot the charge-density-wave amplitude spectral density.
# This channel describes CDW amplitude fluctuations.
ax.plot(
    w_lin.real,
    resolvents.spectral_density(w_lin, "amplitude_CDW"),
    label="CDW",
    ls=linestyles[2],
    c=colors[2]
)

# Plot the longitudinal antiferromagnetic amplitude spectral density.
# The label "l.AFM" stands for longitudinal AFM.
ax.plot(
    w_lin.real,
    resolvents.spectral_density(w_lin, "amplitude_AFM"),
    label="l.AFM",
    ls=linestyles[3],
    c=colors[3]
)

# Plot the transverse antiferromagnetic amplitude spectral density.
# The label "t.AFM" stands for transverse AFM.
ax.plot(
    w_lin.real,
    resolvents.spectral_density(w_lin, "amplitude_AFM_transversal"),
    label="t.AFM",
    ls=linestyles[4],
    c=colors[4]
)

# Shade the continuum region.
resolvents.mark_continuum(ax)

ax.legend()
fig.tight_layout()
plt.show()