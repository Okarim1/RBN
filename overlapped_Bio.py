# -*- coding: utf-8 -*-
"""
@author: okarim
"""
import rbn
import time
import matplotlib.pyplot as plt
import numpy as np
from  tqdm import trange
from scipy import stats

if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    
    start_time = time.time()
    plt.title("O=1")
    plt.ylabel("Fragility")
    plt.xlabel("X")
    
    colors=['b', 'orange', 'g', 'brown', 'purple', 'pink', 'y' ]
    
    red=rbn.RBN()
    
    red.CreateBioNet(1)
    f=[]
    for i in range(100):
        f.append(red.antifragile(100, runs=1, O=1, fraction=1))
        
    plt.errorbar(np.arange(1,red.N+1), np.mean(f,0), label="CD4+ T Cell", 
                 yerr=stats.sem(f, 0), ecolor='r', color=colors[0] )
    
    red.CreateBioNet(2)
    f=[]
    for i in range(100):
        f.append(red.antifragile(100, runs=1, O=1, fraction=1))
    plt.errorbar(np.arange(1,red.N+1), np.mean(f,0), label="Mammalian", 
                 yerr=stats.sem(f, 0), ecolor='r', color=colors[1] )
    
    red.CreateBioNet(3)
    f=[]
    for i in range(100):
        f.append(red.antifragile(100, runs=1, O=1, fraction=1))
    plt.errorbar(np.arange(1,red.N+1), np.mean(f,0), label="Cardiac", 
                 yerr=stats.sem(f, 0), ecolor='r', color=colors[2] )

    red.CreateBioNet(4)
    f=[]
    for i in range(100):
        f.append(red.antifragile(100, runs=1, O=1, fraction=1))
    plt.errorbar(np.arange(1,red.N+1), np.mean(f,0), label="Metabolic", 
                 yerr=stats.sem(f, 0), ecolor='r', color=colors[3] )

    red.CreateBioNet(5)
    f=[]
    for i in range(100):
        f.append(red.antifragile(100, runs=1, O=1, fraction=1))
    plt.errorbar(np.arange(1,red.N+1), np.mean(f,0), label="Death", 
                 yerr=stats.sem(f, 0), ecolor='r', color=colors[4] )
    
    red.CreateBioNet(6)
    f=[]
    for i in range(100):
        f.append(red.antifragile(100, runs=1, O=1, fraction=1))
    plt.errorbar(np.arange(1,red.N+1), np.mean(f,0), label="Arabidopsis", 
                 yerr=stats.sem(f, 0), ecolor='r', color=colors[5] )
    
    red.CreateBioNet(7)
    f=[]
    for i in range(100):
        f.append(red.antifragile(100, runs=1, O=1, fraction=1))
    plt.errorbar(np.arange(1,red.N+1), np.mean(f,0), label="Tumour Cell", 
                 yerr=stats.sem(f, 0), ecolor='r', color=colors[6] )
    
    plt.legend()
    
    print("--- %s seconds ---" % (time.time() - start_time))