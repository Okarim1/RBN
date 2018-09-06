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
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    
    start_time = time.time()
    
    N=100
    p=0.5
    T=100
    
    number_of_iterations=50
    plt.ylabel("fragility")
    plt.xlabel("X")
        
    for K in range(1, 6):
        f=np.zeros(( number_of_iterations, int(N/2) ))
        i=0
        for x in trange(number_of_iterations):
            red=rbn.RBN()
            red.CreateNet(int(K), N, p)
            f[i]=red.antifragile(T, O=1, runs=10)
            i+=1
        g1=np.mean(f, 0)
        plt.plot(np.insert(g1, 0,0), label="K= "+str(K))
        
    plt.legend()
        
    
    print("--- %s seconds ---" % (time.time() - start_time))