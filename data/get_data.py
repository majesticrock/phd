from itertools import groupby
from lib2to3.pytree import convert
import numpy as np
import gzip
import json
import pandas as pd

import os
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

def __to_path(model, subfolder, **kwargs):
    parameters = "/".join(f"{key}={value}" for key, value in kwargs.items())
    return f"/{model}/{subfolder}/{parameters}/"

def __convert_nested_lists_to_arrays__(obj, dtype=None):
    if isinstance(obj, list):
        return np.array([__convert_nested_lists_to_arrays__(sub_obj) for sub_obj in obj], dtype=dtype)
    return obj

def is_ragged(nested_list):
    # Check if all elements in the current list are lists themselves
    if not all(isinstance(i, list) for i in nested_list):
        return False
    first_length = len(nested_list[0])
    
    for sublist in nested_list:
        if not isinstance(sublist, list) or len(sublist) != first_length:
            return True
        if is_ragged(sublist):
            return True
    return False

def convert_nested_lists_to_arrays(obj):
    if isinstance(obj, list):
        dtype = object if is_ragged(obj) else None
        return np.array([__convert_nested_lists_to_arrays__(sub_obj, dtype) for sub_obj in obj], dtype=dtype)
            
    return obj

def lists_to_arrays_panda(panda_obj):
    return panda_obj.applymap(lambda x: convert_nested_lists_to_arrays(x))

def load_panda(model, subfolder, file, **kwargs):
    with gzip.open(f"{CURRENT_DIR}{__to_path(model, subfolder, **kwargs)}{file}", 'rt') as f_open:
        jData = json.load(f_open)
        data = pd.json_normalize(jData, max_level=1)
    return lists_to_arrays_panda(data)

def continuum_params(T, coulomb_scaling, E_F, g, omega_D):
    return {"T": T, "coulomb_scaling": coulomb_scaling, "E_F": E_F, "g": g, "omega_D": omega_D}

def hubbard_params(T, U, V):
    return {"T": T, "U": U, "V": V}

#pd_data = load_panda("modes/square", "test", "resolvents.json.gz", **hubbard_params(0.0, -2.5, -0.1))
#test = pd_data["resolvents.amplitude_AFM"]
#for d, df in pd_data.groupby("Discretization"):
#    print(d)
#    print(df["resolvents.amplitude_AFM"])
#test2 = load_panda("continuum", "test", "gap.json.gz", **continuum_params(0.0, 1.0, 9.3, 12.71, 0.01))
#print(test2)