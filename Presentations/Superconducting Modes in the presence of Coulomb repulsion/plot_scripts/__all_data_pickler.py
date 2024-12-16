import __path_appender as __ap
__ap.append()
import get_data
import pandas as pd
from timeit import timeit

DATA_CUTS = [0, 4, 12]
PICKLE_NAME = "modes/all_data.pkl"
__BASE_FOLDER_TYPE__ = "offset_"
__SC_ONLY_FOLDER__ = "sc_channel_"
__TEST_LOADTIME__ = False

__USE_ONLY_SC_CHANNEL = False

def __print_g_cuts__(dfs, lambda_screening):
    for df in dfs:
        print(df.loc[df["lambda_screening"] == lambda_screening, "g"].max())

def __load_raw__():
    data_5 = get_data.load_all(f"continuum/{__BASE_FOLDER_TYPE__}5/N_k=20000/T=0.0", "resolvents.json.gz").query(
        f"k_F == 4.25 & Delta_max < {DATA_CUTS[len(DATA_CUTS) - 3]}"
        )
    data_10 = get_data.load_all(f"continuum/{__BASE_FOLDER_TYPE__}10/N_k=20000/T=0.0", "resolvents.json.gz").query(
        f"k_F == 4.25 & Delta_max >= {DATA_CUTS[len(DATA_CUTS) - 3]} & Delta_max < {DATA_CUTS[len(DATA_CUTS) - 2]}"
        )
    data_20 = get_data.load_all(f"continuum/{__BASE_FOLDER_TYPE__}20/N_k=20000/T=0.0", "resolvents.json.gz").query(
        f"k_F == 4.25 & Delta_max >= {DATA_CUTS[len(DATA_CUTS) - 2]} & Delta_max < {DATA_CUTS[len(DATA_CUTS) - 1]}"
        )
    data_25 = get_data.load_all(f"continuum/{__BASE_FOLDER_TYPE__}25/N_k=30000/T=0.0", "resolvents.json.gz").query(
        f"k_F == 4.25 & Delta_max >= {DATA_CUTS[len(DATA_CUTS) - 1]}"
        )
    all_data = pd.concat([data_5, data_10, data_20, data_25])
    return all_data

def __load_full_pickle__():
    data_5 = get_data.load_pickle(f"continuum/{__BASE_FOLDER_TYPE__}5", "resolvents.pkl").query(
        f"discretization == 20000 & T == 0.0 & k_F == 4.25 & Delta_max < {DATA_CUTS[len(DATA_CUTS) - 3]}"
        )
    data_10 = get_data.load_pickle(f"continuum/{__BASE_FOLDER_TYPE__}10", "resolvents.pkl").query(
        f"discretization == 20000 & T == 0.0 & k_F == 4.25 & Delta_max >= {DATA_CUTS[len(DATA_CUTS) - 3]} & Delta_max < {DATA_CUTS[len(DATA_CUTS) - 2]}"
        )
    data_20 = get_data.load_pickle(f"continuum/{__BASE_FOLDER_TYPE__}20", "resolvents.pkl").query(
        f"discretization == 20000 & T == 0.0 & k_F == 4.25 & Delta_max >= {DATA_CUTS[len(DATA_CUTS) - 2]} & Delta_max < {DATA_CUTS[len(DATA_CUTS) - 1]}"
        )
    data_25 = get_data.load_pickle(f"continuum/{__BASE_FOLDER_TYPE__}25", "resolvents.pkl").query(
        f"discretization == 30000 & T == 0.0 & k_F == 4.25 & Delta_max >= {DATA_CUTS[len(DATA_CUTS) - 1]}"
        )
    
    __arr__ = [data_5, data_10, data_20]
    print("1. cuts: ")
    __print_g_cuts__(__arr__, 0)
    print("2. cuts: ")
    __print_g_cuts__(__arr__, 1)
    print("3. cuts: ")
    __print_g_cuts__(__arr__, 1e-4)

    if __USE_ONLY_SC_CHANNEL:
        sc_data_10 = get_data.load_pickle(f"continuum/{__SC_ONLY_FOLDER__}10", "resolvents.pkl").query(
            f"discretization == 20000 & T == 0.0 & k_F == 4.25 & Delta_max >= {DATA_CUTS[len(DATA_CUTS) - 3]} & Delta_max < {DATA_CUTS[len(DATA_CUTS) - 2]}"
            )
        sc_data_20 = get_data.load_pickle(f"continuum/{__SC_ONLY_FOLDER__}20", "resolvents.pkl").query(
            f"discretization == 20000 & T == 0.0 & k_F == 4.25 & Delta_max >= {DATA_CUTS[len(DATA_CUTS) - 2]} & Delta_max < {DATA_CUTS[len(DATA_CUTS) - 1]}"
            )
        sc_data_25 = get_data.load_pickle(f"continuum/{__SC_ONLY_FOLDER__}25", "resolvents.pkl").query(
            f"discretization == 30000 & T == 0.0 & k_F == 4.25 & Delta_max >= {DATA_CUTS[len(DATA_CUTS) - 1]}"
            )
        
        data_5.query("coulomb_scaling == 0", inplace=True)
        data_10.query("coulomb_scaling == 0", inplace=True)
        data_20.query("coulomb_scaling == 0", inplace=True)
        data_25.query("coulomb_scaling == 0", inplace=True)
        
        return pd.concat([ data_5, data_10, data_20, data_25, sc_data_10, sc_data_20, sc_data_25])
        

    all_data = pd.concat([ data_5, data_10, data_20, data_25])
    return all_data

def load_and_pickle():
    all_data = __load_full_pickle__()
    #all_data = __load_raw__()
    all_data.to_pickle(PICKLE_NAME)
    return all_data
    
def load_pickled():
    return pd.read_pickle(PICKLE_NAME)

if __name__ == "__main__":
    load_and_pickle()
    if __TEST_LOADTIME__:        
        print("Loading jsons... ")
        print(timeit(__load_raw__, number=10))
        print("Loading pickle... ")
        print(timeit(load_pickled, number=10))
        print("Loading large pickle... ")
        print(timeit(__load_full_pickle__, number=10))