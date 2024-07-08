import matplotlib.pyplot as plt
import numpy as np
import __path_appender as __ap
__ap.append()
from get_data import *

MEV_FACTOR = 1e3

fig, ax = plt.subplots()

pd_data = load_panda("continuum", "test", "gap.json.gz", 
                     **continuum_params(0.0, 1.0, 9.3, 10., 0.01)).iloc[0]
pd_data["ks"] -= pd_data["k_F"]

ax.plot(pd_data["ks"], label="Discretization")

ax.set_xlabel(r"$n$")
ax.set_ylabel(r"$k / k_\mathrm{F}$")
ax.legend()
fig.tight_layout()

plt.show()