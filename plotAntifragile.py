# -*- coding: utf-8 -*-
"""
@author: okarim
"""
import rbn
import time
import matplotlib.pyplot as plt
import numpy as np
from  tqdm import trange

if __name__ == '__main__':
    """
    Plots the antifragility of RBNs varying perturbations
    """
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    start_time = time.time()
    
    K=2.0
    N=100
    p=0.5
    T=100

    number_of_iterations=50
    FR= np.zeros((5, int(N/2)))
    for K in range(1, 6):
        f=np.zeros(( number_of_iterations, int(N/2) ))
        i=0
        for x in trange(number_of_iterations):
            red=rbn.RBN()
            red.CreateNet(int(K), N, p)
            f[i]=red.antifragile(T, O=1, runs=10)
            i+=1
        plt.figure()
        g1=np.mean(f, 0)
        plt.plot(np.insert(g1, 0,0))
        plt.title("K= "+str(K))
        plt.ylabel("Fragility")
        plt.xlabel("X")
        plt.show()
        
        plt.figure()
        plt.boxplot(f)
        plt.title("K= "+str(K))
        plt.ylabel("Fragility")
        plt.xlabel("X")
        plt.show()
        
        FR[K-1]=g1
    
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#    ax.plot_surface(FR[:, 0], FR[:, 1], )
#    plt.show()
