import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

data_folder = "data/T0.5/"

CDW = (np.loadtxt(data_folder + "cdw.txt")).transpose()
SC =  (np.loadtxt(data_folder + "sc.txt")).transpose()
ETA = (np.loadtxt(data_folder + "eta.txt")).transpose()

labels = ["T", "U"]
T_SIZE = len(CDW)
U_SIZE = len(CDW[0])

with open(data_folder + "cdw.txt") as fp:
    for i, line in enumerate(fp):
        if i == 2:
            ls = line.split()
            labels[0] = ls[1].split("_")[0]
            U = np.linspace(float(ls[1].split("=")[1]), float(ls[2].split("=")[1]), U_SIZE+1)[:U_SIZE]
        elif i == 3:
            ls = line.split()
            labels[1] = ls[1].split("_")[0]
            T = np.linspace(float(ls[1].split("=")[1]), float(ls[2].split("=")[1]), T_SIZE+1)[:T_SIZE]
        elif i > 3:
            break

def func(T, a, b):
    ret_arr = np.zeros(len(T))
    for i in range(0, len(T)):
        if(b - T[i] > 0):
            ret_arr[i] = a * np.sqrt(b - T[i])
        else:
            ret_arr[i] = 0
    return ret_arr
    

fitted_cdw = np.zeros(U_SIZE)
fitted_sc = np.zeros(U_SIZE)

#CDW = np.sqrt(CDW*CDW + ETA*ETA + SC*SC)

for i in range(0, U_SIZE):
    RANGE = 10
    END = 0
    guess = np.zeros(2)
    skip = True
    for j in range(0, T_SIZE):
        if(CDW[i][j] > 1e-8):
            skip = False
            if(guess[0] == 0):
                guess[0] = CDW[i][j]
        elif(not skip):
            guess[1] = T[j]
            END = j + 1
            if END < 2:
                skip = True
            break

    if(not skip and END != 0):
        while(RANGE > END):
            RANGE -= 1
        popt, pcov =  curve_fit(func, T[END-RANGE:END], CDW[i][END-RANGE:END], guess)
        fitted_cdw[i] = popt[1]

    RANGE = 10
    END = 0
    guess = np.zeros(2)
    skip = True
    for j in range(0, T_SIZE):
        if(SC[i][j] > 1e-8):
            skip = False
            if(guess[0] == 0):
                guess[0] = SC[i][j]
        elif(not skip):
            guess[1] = T[j]
            END = j + 1
            if END < 2:
                skip = True
            break
    
    if(not skip and END != 0):
        while(RANGE > END):
            RANGE -= 1
        popt, pcov =  curve_fit(func, T[END-RANGE:END], SC[i][END-RANGE:END], guess)
        fitted_sc[i] = popt[1]


plt.figure()
plt.plot(U, fitted_sc, "k")

plt.fill_between(U, fitted_sc, 0, alpha=0.4, label="SC")
plt.plot(U, fitted_cdw, "k")
plt.fill_between(U, fitted_cdw, 0, alpha=0.4, label="CDW")
plt.xlabel("$" + labels[0] + " / t$")
plt.ylabel("$T / t$")
plt.legend(loc="upper center")

import os
plt.savefig(f"python/build/{os.path.basename(__file__).split('.')[0]}.pdf")
plt.show()