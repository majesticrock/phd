import os, sys

def append():
    __base = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # parent folder of "phd"
    sys.path.append(__base + "/PhdUtility/python")