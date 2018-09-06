# -*- coding: utf-8 -*-
"""
@author: omark
"""
import rbn
import time
import matplotlib.pyplot as plt
import numpy as np
from  tqdm import trange

if __name__ == '__main__':
    """
     Plots the minimum fragility for each % of attractors against perturbations and fragility value
    """
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    start_time = time.time()
    
    N=5
    p=0.5
    K=2.0
    T=100
    
    runs=1000
    
    vf = []
    vp = []
    numAtt=[]
    avLen=[]
    
    red=rbn.RBN()
    
    for i in trange(runs):
        red.CreateNet(K, N, p)
        f=red.antifragile(T, O=1, runs=100)
        vf.append(np.amin(f))
        vp.append(np.argmin(f))
        A = red.Attractors(T, runs=0)
        numAtt.append(len(A))
        if(len(A) == 0):
            avLen.append(T)
        else:
            edos=0
            for j in A:
                edos+=len(j)
            avLen.append(edos/len(A))

    plt.plot(numAtt, vf, 'o')
    plt.ylabel("Fragility")
    plt.xlabel("Attractors")
    plt.show()  
    
    plt.plot(numAtt, vp, 'o')
    plt.ylabel("Perturbations")
    plt.xlabel("Attractors")
    plt.show()
    
    plt.plot(avLen, vf, 'o')
    plt.ylabel("Fragility")
    plt.xlabel("Attractors Length")
    plt.show()  
    
    plt.plot(avLen, vp, 'o')
    plt.ylabel("Perturbations")
    plt.xlabel("Attractors Length")
    plt.show()
    