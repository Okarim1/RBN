# -*- coding: utf-8 -*-
"""
@author: okarim
"""
import rbn
import time
import numpy as np
from  tqdm import trange

def fragility(red, X, O, T, runs):
    """
    P=Probability of generating antifragile networks
    AF=average fragility
    AIC=average initial complexity
    AFC=average final complexity
    """
    
    AF=0
    AIC=0
    AFC=0
    P=0
    for i in range(runs):
        initial = np.random.randint(0, 2, red.N)
        State=red.RunNet(T, initial)
        C0 = rbn.complexity(State)
        
        State=red.RunNet(T, initial, X, O)
        C = rbn.complexity(State)
        
        f=rbn.fragility(C, C0, X, O, red.N, T)
        
        AF+=f
        AIC+=np.mean(C0)
        AFC+=np.mean(C)
        if f < 0:
            P+=1
        
    return AF/runs, AIC/runs, AFC/runs, P/runs

if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    
    start_time = time.time()
    
    N=100
    p=0.5
    T=1000
    
    number_of_iterations=100
        
    for K in range(1, 5):
        AF=0
        AIC=0
        AFC=0
        AP=0
        for x in trange(number_of_iterations):
            red=rbn.RBN()
            red.CreateNet(int(K), N, p)
            f, C0, C, P =fragility(red, X=20, O=1, T=1000, runs=10)
            
            AF+=f
            AIC+=C0
            AFC+=C
            AP+=P
        
        AF/=number_of_iterations
        AIC/=number_of_iterations
        AFC/=number_of_iterations
        AP/=number_of_iterations
        
        print("K="+str(K)+":")
        print("Probability of generating antifragile networks= "+str(AP))
        print("Average fragility= "+str(AF))
        print("Average initial complexity= "+str(AIC))
        print("Average final complexity= "+str(AFC)+"\n")
    
    print("--- %s seconds ---" % (time.time() - start_time))