# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 17:41:13 2018

@author: omark
"""
import rbn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def CreateBioNet():
    red=rbn.RBN()
    red.N=18
    
    data=pd.read_csv("BioConInv.csv", sep=",", header=-1)
    data=data.fillna(0)
    x=np.array(data.loc[:, 1:]).astype(int)
    red.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x])
    
    data=pd.read_csv("BioBool.csv", sep=",", header=0)
    data=data.fillna(0)
    y=np.array(data.loc[:, :]).astype(int).transpose()
    red.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y])
    
    red.K=np.count_nonzero(red.Con)/12
    return red
    
if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    red=rbn.RBN()
    red.CreateBioNet()
    
    State=red.RunNet(20, Bio=True)
    plt.imshow(State, cmap='Greys', interpolation='None')
    
    f=red.antifragile(100, runs=500, O=1, fraction=1)
    
    plt.show()
    plt.plot(f, label="BioNet")
    plt.title("O=1")
    plt.ylabel("Fragility")
    plt.xlabel("X")
    
    p=150/(3242+150)
    r=50
    f=np.zeros(( r,  18))
    for i in range(r):
        red.CreateNet(6.5, 18, p)
        f[i]=red.antifragile(100, runs=10, O=1, fraction=1)
    
    plt.plot(np.mean(f, 0), label="RBN")
    plt.title("O=1")
    plt.ylabel("Fragility")
    plt.xlabel("X")
    
    plt.legend()
    