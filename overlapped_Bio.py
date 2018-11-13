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
    plt.title("O=1")
    plt.ylabel("Fragility")
    plt.xlabel("X")
    
    red=rbn.RBN()
    
    red.CreateBioNet(1)
    f=red.antifragile(100, runs=500, O=1, fraction=1)
    plt.plot(np.arange(1,red.N+1), f, label="CD4+ T Cell")
    
    red.CreateBioNet(2)
    f=red.antifragile(100, runs=500, O=1, fraction=1)
    plt.plot(np.arange(1,red.N+1), f, label="Mammalian")
    
    red.CreateBioNet(3)
    f=red.antifragile(100, runs=500, O=1, fraction=1)
    plt.plot(np.arange(1,red.N+1), f, label="Cardiac")

    red.CreateBioNet(4)
    f=red.antifragile(100, runs=500, O=1, fraction=1)
    plt.plot(np.arange(1,red.N+1), f, label="Metabolic")

    red.CreateBioNet(5)
    f=red.antifragile(100, runs=500, O=1, fraction=1)
    plt.plot(np.arange(1,red.N+1), f, label="Death")
    
    red.CreateBioNet(6)
    f=red.antifragile(100, runs=500, O=1, fraction=1)
    plt.plot(np.arange(1,red.N+1), f, label="Arabidopsis")
    
    red.CreateBioNet(7)
    f=red.antifragile(100, runs=500, O=1, fraction=1)
    plt.plot(np.arange(1,red.N+1), f, label="Tumour Cell")
    
    plt.legend()
    
    print("--- %s seconds ---" % (time.time() - start_time))