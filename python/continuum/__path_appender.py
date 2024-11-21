import os, sys

def __base():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) #  "phd"

def append():
    sys.path.append(os.path.join(__base(), "PhdUtility", "python"))
    sys.path.append(os.path.join(os.path.dirname(__base()), "raw_data_phd"))
    sys.path.append(os.path.join(__base(), "data"))
    