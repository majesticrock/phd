import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Linearly fits to a data set and plots the fit, returns the best fit parameters
def ez_linear_fit(x_data, y_data, plotter=plt, x_space=None, **plotter_kwargs):
    def func(x, a, b):
        return a*x + b
    popt, pcov = curve_fit(func, x_data, y_data)
    
    if not hasattr(x_space, "__len__"):
        x_space = np.linspace(min(x_data), max(x_data), 200)
    plotter.plot(x_space, func(x_space, *popt), **plotter_kwargs)
    
    return popt, pcov