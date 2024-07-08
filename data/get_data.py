import numpy as np
import gzip
import json
import pandas as pd

import os
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

def __to_path(model, subfolder, **kwargs):
    parameters = "/".join(f"{key}={value}" for key, value in kwargs.items())
    return f"/{model}/{subfolder}/{parameters}/"

def convert_lists_to_arrays(obj):
    if isinstance(obj, list):
        return np.array(obj)
    if isinstance(obj, dict):
        return {k: convert_lists_to_arrays(v) for k, v in obj.items()}
    return obj

def __load_panda__(file):
    with gzip.open(file, 'rt') as f_open:
        jData = json.load(f_open, object_hook=convert_lists_to_arrays)
        data = pd.json_normalize(jData, max_level=1)
    return data

def load_panda(model, subfolder, file, **kwargs):
    return __load_panda__(f"{CURRENT_DIR}{__to_path(model, subfolder, **kwargs)}{file}")

def continuum_params(T, coulomb_scaling, E_F, g, omega_D):
    return {"T": T, "coulomb_scaling": coulomb_scaling, "E_F": E_F, "g": g, "omega_D": omega_D}

def hubbard_params(T, U, V):
    return {"T": T, "U": U, "V": V}

import glob
def load_all(folder, pattern):
    all_files = glob.glob(os.path.join(f"{CURRENT_DIR}/{folder}", f"**/{pattern}"), recursive=True)
    return pd.concat((__load_panda__(file) for file in all_files), ignore_index=True)

#pd_data = load_all("continuum/test", "gap.json.gz")
#print(pd_data)

#pd_data = load_panda("modes/square", "test", "resolvents.json.gz", **hubbard_params(0.0, -2.5, -0.1))
#test = pd_data["resolvents.amplitude_AFM"]
#for d, df in pd_data.groupby("Discretization"):
#    print(d)
#    print(df["resolvents.amplitude_AFM"])
#test2 = load_panda("continuum", "test", "gap.json.gz", **continuum_params(0.0, 1.0, 9.3, 12.71, 0.01))
#print(test2)