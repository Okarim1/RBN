# -*- coding: utf-8 -*-
"""
@author: okarim
"""
import rbn
import time
import matplotlib.pyplot as plt
import numpy as np
from  tqdm import tqdm

if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    
    start_time = time.time()
    
    N=100
    p=0.5
    T=100
    
    number_of_iterations=1000
    plt.ylabel("Initial Complexity")
    plt.xlabel("K")
    red=rbn.RBN()
    g1=[]
    
    rango=np.arange(0.1, 5.0, 0.1)
    for K in tqdm(rango):
        C=np.zeros(number_of_iterations)
        i=0
        for x in range(number_of_iterations):
            red.CreateNet(K, N, p)
            State = red.RunNet(2*T)
            C[i]=np.mean(rbn.complexity(State[-T:]))
            i+=1
        g1.append(np.mean(C))
    
    plt.plot(np.arange(0.1,5.0,0.1), g1, label="K= "+str(K))
    plt.title("Average Complexity")    
    
    print("--- %s seconds ---" % (time.time() - start_time))