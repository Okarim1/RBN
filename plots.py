# -*- coding: utf-8 -*-
"""
@author: okarim
"""
import rbn
import time
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import hamming

if __name__ == '__main__':
    """
    Plot examples of a network with and without perturbations
    """
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    start_time = time.time()
    
    K=2.0
    N=20
    p=0.5
    T=40
    
    X=2 # how many perturbations
    O=1 # how often the perturbations take place
    
    red=rbn.RBN()
    red.CreateNet(K, N, p)
    
    initial = np.random.randint(0, 2, red.N)
    State=red.RunNet(T, initial)
    plt.imshow(State, cmap='Greys', interpolation='None')
    #, aspect='auto'
    plt.show()
    C1 = rbn.complexity(State)
    print("Complejidad: "+str(np.mean(C1)))
    initial2 = np.random.randint(0, 2, red.N)
    State2=red.RunNet(T, initial2)
    plt.imshow(State2, cmap='Greys', interpolation='None')
    plt.show()
    C3 = rbn.complexity(State2)
    print("Complejidad: "+str(np.mean(C3)))
    print("\nDistancia inicial: ")
    print(hamming(State[0], State2[0])) #distance
    print("Distancia final: ")
    print(hamming(State[T-1], State2[T-1])) #distance
    
    State=red.RunNet(T, initial, X, O)
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.show()
    C2 = rbn.complexity(State)
    print("Complejidad: "+str(np.mean(C2)))
    print("Fragilidad: "+str(rbn.fragility(C2, C1, X, O, red.N, T))) 
    State2=red.RunNet(T, initial2, X, O)
    plt.imshow(State2, cmap='Greys', interpolation='None')
    plt.show()
    C4 = rbn.complexity(State2)
    print("Complejidad: "+str(np.mean(C4)))
    print("Fragilidad: "+str(rbn.fragility(C4, C3, X, O, red.N, T)))
    print("\nDistancia inicial: ")
    print(hamming(State[0], State2[0])) #distance
    print("Distancia final: ")
    print(hamming(State[T-1], State2[T-1])) #distance
    
    f=red.antifragile(T, O=1, runs=100)
    plt.plot(f)
    plt.ylabel("Fragility")
    plt.xlabel("Perturbations")
    plt.show()  
    