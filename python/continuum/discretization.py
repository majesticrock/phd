import matplotlib.pyplot as plt
import numpy as np
import __path_appender as __ap
__ap.append()
from get_data import *

MEV_FACTOR = 1e3

fig, ax = plt.subplots()

main_df = load_panda("continuum", "test", "gap.json.gz", 
                     **continuum_params(0.0, 0.0, 9.3, 5., 10.))
pd_data = main_df["data"]
pd_data["ks"] /= main_df["k_F"]

ax.plot(pd_data["ks"], label="Discretization")

ax.set_xlabel(r"$n$")
ax.set_ylabel(r"$k / k_\mathrm{F}$")
ax.legend()
fig.tight_layout()
plt.grid()

plt.show()