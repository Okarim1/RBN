# -*- coding: utf-8 -*-
"""
@author: okarim
"""
import rbn
import time
import matplotlib.pyplot as plt
import numpy as np
from  tqdm import tqdm
from scipy import stats

if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    
    start_time = time.time()
    
    N=100
    p=0.5
    T=1000
    
    number_of_iterations=1000
    plt.ylabel("Initial Complexity")
    plt.xlabel("K")
    red=rbn.RBN()
    g1=[]
    yerr=[]
    
    rango=np.arange(0.1, 5.1, 0.1)
    for K in tqdm(rango):
        C=[]
        for x in range(number_of_iterations):
            red.CreateNet(K, N, p)
            State = red.RunNet(2*T)
            C.append(np.mean(rbn.complexity(State[-T:])))
        g1.append(np.mean(C))
        yerr.append(C)
    
    plt.errorbar(rango, g1, label="K= "+str(K), yerr=stats.sem(yerr,1), ecolor='r')
    plt.title("Average Complexity")    
    plt.ylim(top=1)
    plt.savefig("Figure_4a.eps")
    np.savez("initComplex.npz", g1)
    print("--- %s seconds ---" % (time.time() - start_time))