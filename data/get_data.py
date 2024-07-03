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
    return panda_obj.apply(lambda x: convert_nested_lists_to_arrays(x))

def dicts_to_panda(panda_obj):
    return panda_obj.apply(lambda x: pd.DataFrame.from_dict(x) if isinstance(x, dict) else x)

def load_panda_series(model, subfolder, file, **kwargs):
    with open(f"{CURRENT_DIR}{__to_path(model, subfolder, **kwargs)}{file}", 'rt') as f_open:
        #data = pd.concat([pd.read_json(f_open, typ="series"), pd.Series(data=kwargs)])
        #data = pd.read_json(f_open, typ="series").groupby("resolvents")
        jData = json.load(f_open)
        data = pd.json_normalize(jData, max_level=1)
        data2 = pd.json_normalize(jData, max_level=1)
        data2.at[0, "Discretization"] = 1000
        
    return pd.concat([data, data2]).reset_index()
    return lists_to_arrays_panda(data)

def continuum_params(T, E_F, g, omega_D):
    return {"T": T, "E_F": E_F, "g": g, "omega_D": omega_D}

def hubbard_params(T, U, V):
    return {"T": T, "U": U, "V": V}

pd_data = load_panda_series("modes/square", "test", "resolvents.json", **hubbard_params(0.0, -2.5, -0.1))
test = pd_data["resolvents.amplitude_AFM"]
for d, df in pd_data.groupby("Discretization"):
    print(d)
    print( df["resolvents.amplitude_AFM"])
# for x in test:
#     print(x)

#print(load_panda_series("continuum", "test", "gap.json.gz", **{"T": 0.0, "E_F": 9.3, "g": 12.71, "omega_D": 0.01}))