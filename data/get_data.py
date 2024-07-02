import numpy as np
import gzip
import pandas as pd


with gzip.open(f"TODO.json.gz", 'rt') as f_open:
    data = pd.read_json(f_open)