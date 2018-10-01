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
    Plots the minimum fragility for each K against perturbations and fragility value
    """
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    start_time = time.time()
    
    N=100
    p=0.5
    T=100
    runs=50
    
    red=rbn.RBN()
    
    vfm = []
    vpm = []
    
    vfbp = []
    vpbp = []
    
    rango=np.arange(0.1, 5.0, 0.1)
    for K in tqdm(rango):
        vf = []
        vp = []
        for i in range(runs):
            red.CreateNet(K, N, p)
            f=red.antifragile(T, X=40, runs=10)
            vf.append(np.amin(f))
            vp.append(np.argmin(f))
        vfm.append(np.mean(vf))
        vpm.append(np.mean(vp))
        
        vfbp.append(vf)
        vpbp.append(vp)
        
    plt.plot(rango, vfm)
    plt.ylabel("Fragility")
    plt.xlabel("K")
    plt.show()   
    plt.plot(rango, vpm)
    plt.ylabel("O")
    plt.xlabel("K")
    plt.show() 
    
    plt.boxplot(vfbp)
    plt.ylabel("Fragility")
    plt.xlabel("K")
    plt.show()
    plt.boxplot(vpbp)
    plt.ylabel("O")
    plt.xlabel("K")
    plt.show() 
    
    np.savez("va.npz", vfm, vpm, vfbp, vpbp)
    
    #x=np.arange(0.1, 3.0, 0.1)
    #y=vfm
    #print(np.polyfit(np.log(x), y, 1))
    #print(np.polyfit(x, np.log(y), 1, w=np.sqrt(y)))