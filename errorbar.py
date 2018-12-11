# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 18:56:22 2018

@author: Okarim
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
 
data=np.load("npz/newX.npz")
vfm=data['arr_0']
vpm=data['arr_1']
vfbp=data['arr_2']
vpbp =data['arr_3']

rango=np.arange(0.1, 5.0, 0.1)

yerr=stats.sem(vfbp,1)
plt.errorbar(rango, vfm, yerr=yerr, ecolor='r')
plt.ylabel("Fragility")
plt.xlabel("K")
plt.show()

yerr=stats.sem(vpbp,1)
plt.errorbar(rango, vpm, yerr=yerr, ecolor='r')
plt.ylabel("X")
plt.xlabel("K")
plt.show() 

