import os, sys

def __base():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # parent folder of "phd"

def append():
    sys.path.append(__base() + "/PhdUtility/python")
    sys.path.append(__base() + "/phd/data")