# -*- coding: utf-8 -*-
"""
@author: omark
"""
import rbn
import time
import matplotlib.pyplot as plt
import numpy as np
from  tqdm import tqdm

if __name__ == '__main__':
    """
    Plots the minimum fragility for each M against perturbations and fragility value
    """
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    start_time = time.time()
    
    K=2
    N=100
    p=0.5
    T=100
    
    Prob=0.9

    runs=50
    red=rbn.RBN()
    
    vfm = []
    vpm = []
    
    vfbp = []
    vpbp = []
    
    rango=np.arange(-10, 10)
    for M in tqdm(rango):
        vf = []
        vp = []
        for i in range(runs):
            if M in [-1, 0, 1]:
                red.CreateNet(K, N, p)
            else:
                red.CreateNetMod(K, N, p, M, Prob)
            f=red.antifragile(T, O=1, runs=100)
            vf.append(np.amin(f))
            vp.append(np.argmin(f))
        vfm.append(np.mean(vf))
        vpm.append(np.mean(vp))
        
        vfbp.append(vf)
        vpbp.append(vp)
        
    plt.plot(rango, vfm)
    plt.ylabel("Fragility")
    plt.xlabel("M")
    plt.show()   
    plt.plot(rango, vpm)
    plt.ylabel("Perturbations")
    plt.xlabel("M")
    plt.show() 
    
    plt.boxplot(vfbp)
    plt.ylabel("Fragility")
    plt.xlabel("M")
    plt.show()
    plt.boxplot(vpbp)
    plt.ylabel("Perturbations")
    plt.xlabel("M")
    plt.show() 