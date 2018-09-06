# -*- coding: utf-8 -*-
"""
@author: okarim
"""
import rbn
import time
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from  tqdm import trange

if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    
    start_time = time.time()
    
    N=100
    p=0.5
    T=100
    K=1
    
    maxO=50
    number_of_iterations=50
    Z= np.zeros([maxO,50])
    j=0
    for O in trange(maxO):
        f=np.zeros(( number_of_iterations, int(N/2) ))
        i=0
        for x in range(number_of_iterations):
            red=rbn.RBN()
            red.CreateNet(K, N, p)
            f[i]=red.antifragile(T, runs=10, O=O+1)
            i+=1
            
        Z[j,:]=np.mean(f, 0)
        j+=1
    
    
    
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    X = np.arange(1, 51)
    Y = np.arange(1, maxO+1)
    X, Y = np.meshgrid(X, Y)
    
    
    ax.plot_surface(X, Y, Z)
    cset = ax.contour(X, Y, Z, zdir='z', offset=-0.2, cmap=cm.coolwarm)
    cset = ax.contour(X, Y, Z, zdir='x', offset=0, cmap=cm.coolwarm)
    cset = ax.contour(X, Y, Z, zdir='y', offset=50, cmap=cm.coolwarm)
    
    ax.set_xlabel('X')
    ax.set_xlim(0, 50)
    ax.set_ylabel('O')
    ax.set_ylim(0, 50)
    ax.set_zlabel('fragility')
    ax.set_zlim(-0.2, 0.2)
    
    plt.show()
