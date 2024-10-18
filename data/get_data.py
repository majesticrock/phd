import numpy as np
import gzip
import json
import pandas as pd

import os
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

def __to_path(model, subfolder, **kwargs):
    params = os.path.join(*(f"{key}={value}" for key, value in kwargs.items()))
    return os.path.join(CURRENT_DIR, model, subfolder, params)

def convert_lists_to_arrays(obj):
    if isinstance(obj, list):
        return np.array(obj)
    if isinstance(obj, dict):
        return {k: convert_lists_to_arrays(v) for k, v in obj.items()}
    return obj

def __load_panda__(file):
    with gzip.open(file, 'rt') as f_open:
        jData = json.load(f_open, object_hook=convert_lists_to_arrays)
    
    if 'data' in jData:
        main_data = {k: v for k, v in jData.items() if k != 'data'}
        main_df = pd.DataFrame([main_data])
        main_df['data'] = [pd.DataFrame(jData['data'])]
    else:
        main_df = pd.json_normalize(jData, max_level=1)
    return main_df

def load_panda(model, subfolder, file, **kwargs):
    data = __load_panda__(os.path.join(__to_path(model, subfolder, **kwargs), file)).iloc[0]
    print(f"Loaded data has been produced on {data['time']}")
    return data

def continuum_params(N_k, T, coulomb_scaling, screening, k_F, g, omega_D):
    cp_screening = 0.0 if abs(coulomb_scaling) < 1e-12 else screening
    return {"N_k" : N_k , "T": T, "coulomb_scaling": coulomb_scaling, "screening": cp_screening, "k_F": k_F, "g": g, "omega_D": omega_D}

def hubbard_params(T, U, V):
    return {"T": T, "U": U, "V": V}

def __all_files__(folder, pattern):
    if isinstance(folder, str):
        folder = os.path.normpath(folder)
    search_obj = os.path.join(CURRENT_DIR, folder, os.path.normpath("**/"), pattern)
    return glob.glob(search_obj, recursive=True)

import glob
def load_all(folder, pattern, condition=None):
    all_files = __all_files__(folder, pattern)
    if condition is not None:
        if hasattr(condition, "__len__"):
            for cond in condition:
                if os.name == 'nt':
                    cond = f"\\{cond}\\"
                else:
                    cond = f"/{cond}/"
                all_files = [file for file in all_files if cond in file]
        else:
            if os.name == 'nt':
                cond = f"\\{condition}\\"
            else:
                cond = f"/{condition}/"
            all_files = [file for file in all_files if condition in file]
    return pd.concat((__load_panda__(file) for file in all_files), ignore_index=True)

#pd_data = load_all("continuum/disc_2000", "gap.json.gz")
#print(pd_data["data"][0]["Delta_Coulomb"] + pd_data["data"][0]["Delta_Phonon"])