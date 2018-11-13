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
    fraction=1
    plt.ylabel("Final Complexity")
    plt.xlabel("X")
    red=rbn.RBN()
    p1=np.zeros((5, int(N/fraction)))
    
    for K in range(1, 6):
        f=np.zeros(( number_of_iterations, int(N/fraction) ))
        i=0
        for x in trange(number_of_iterations):
            red.CreateNet(int(K), N, p)
            f[i]=red.antifragile(T, O=1, runs=10, fraction=fraction)
            i+=1
        g1=np.mean(f, 0)
        p1[K-1]=np.sum(np.array(f) < 0, axis=0)/number_of_iterations
        #g1=np.insert(g1, 0,0)
        plt.plot(np.arange(1,int(N/fraction)+1), g1, label="K= "+str(K))
    plt.title("Average Complexity")
    plt.legend()
    
    plt.show()
    plt.ylabel("Probability")
    for K in range(1, 5):
        plt.plot(np.arange(1,int(N/fraction)+1), p1[K-1], label="K= "+str(K))
    plt.title("Probility of generating antifragile networks")
    plt.legend()
    
#    red.CreateBioNet()
    
#    f=red.antifragile(100, runs=500, O=1, fraction=1)
#    plt.plot(f, label="BioNet")
#    
    
    print("--- %s seconds ---" % (time.time() - start_time))