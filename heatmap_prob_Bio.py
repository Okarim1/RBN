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
    
    red=rbn.RBN()
    red.CreateBioNet(6)
    N=red.N
    T=100
    
    maxO=20
    number_of_iterations=1
    fraction=1
    
    Z= np.zeros([maxO,int(N/fraction)])
    j=0
    for O in trange(maxO):
        f=np.zeros(( number_of_iterations, int(N/fraction) ))
        i=0
        for x in range(number_of_iterations):
            f[i]=red.antifragile(T, runs=500, O=O+1)
            i+=1
            
        Z[j,:]=np.sum(np.array(f), axis=0)/number_of_iterations
        j+=1
    
    
    np.savez("BioProb6.npz", Z)
    
    fig = plt.figure()
    
    fig, ax = plt.subplots()
    im = ax.imshow(Z, extent=[1,red.N,maxO,1])
    plt.xlabel('X')
    plt.ylabel('O')
    plt.title("Probability with biological Boolean network")
    fig.colorbar(im)