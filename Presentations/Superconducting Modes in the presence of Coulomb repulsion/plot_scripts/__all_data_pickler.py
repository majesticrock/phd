import __path_appender as __ap
__ap.append()
from get_data import *
import pandas as pd
from timeit import timeit

DATA_CUTS = [0, 4, 12]
PICKLE_NAME = "modes/all_data.pkl"

__TEST_LOADTIME__ = False

def __load_raw__():
    data_5 = load_all("continuum/offset_5/N_k=20000/T=0.0", "resolvents.json.gz").query(
        f"k_F == 4.25 & Delta_max < {DATA_CUTS[len(DATA_CUTS) - 3]}"
        )
    data_10 = load_all("continuum/offset_10/N_k=20000/T=0.0", "resolvents.json.gz").query(
        f"k_F == 4.25 & Delta_max >= {DATA_CUTS[len(DATA_CUTS) - 3]} & Delta_max < {DATA_CUTS[len(DATA_CUTS) - 2]}"
        )
    data_20 = load_all("continuum/offset_20/N_k=20000/T=0.0", "resolvents.json.gz").query(
        f"k_F == 4.25 & Delta_max >= {DATA_CUTS[len(DATA_CUTS) - 2]} & Delta_max < {DATA_CUTS[len(DATA_CUTS) - 1]}"
        )
    data_25 = load_all("continuum/offset_25/N_k=30000/T=0.0", "resolvents.json.gz").query(
        f"k_F == 4.25 & Delta_max >= {DATA_CUTS[len(DATA_CUTS) - 1]}"
        )
    all_data = pd.concat([data_5, data_10, data_20, data_25])
    return all_data

def load_and_pickle():
    all_data = __load_raw__()
    all_data.to_pickle(PICKLE_NAME)
    
def load_pickled():
    return pd.read_pickle(PICKLE_NAME)

if __name__ == "__main__":
    load_and_pickle()
    if __TEST_LOADTIME__:
        print("Loading jsons... ")
        print(timeit(__load_raw__, number=10))
        print("Loading pickle... ")
        print(timeit(load_pickled, number=10))