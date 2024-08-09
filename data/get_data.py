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
    
    if 'data' in jData:
        main_data = {k: v for k, v in jData.items() if k != 'data'}
        main_df = pd.DataFrame([main_data])
        main_df['data'] = [pd.DataFrame(jData['data'])]
    else:
        main_df = pd.json_normalize(jData, max_level=1)
    return main_df

def load_panda(model, subfolder, file, **kwargs):
    data = __load_panda__(f"{CURRENT_DIR}{__to_path(model, subfolder, **kwargs)}{file}").iloc[0]
    print(f"Loaded data has been produced on {data['time']}")
    return data

def continuum_params(T, coulomb_scaling, k_F, g, omega_D):
    return {"T": T, "coulomb_scaling": coulomb_scaling, "k_F": k_F, "g": g, "omega_D": omega_D}

def hubbard_params(T, U, V):
    return {"T": T, "U": U, "V": V}

import glob
def load_all(folder, pattern):
    all_files = glob.glob(os.path.join(f"{CURRENT_DIR}/{folder}", f"**/{pattern}"), recursive=True)
    return pd.concat((__load_panda__(file) for file in all_files), ignore_index=True)

#pd_data = load_all("continuum/disc_2000", "gap.json.gz")
#print(pd_data["data"][0]["Delta_Coulomb"] + pd_data["data"][0]["Delta_Phonon"])