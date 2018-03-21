# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 21:06:32 2018

@author: Luis
"""

import matplotlib.pyplot as plt
import numpy as np
import json

type(np.int32(0).item())
type(np.asscalar(np.int32(0)))

def conv(lista):
    cadena = json.dumps(lista)
    newstring = cadena.replace("[","")
    new2 = newstring.replace("]","")
    new3 = new2.replace(",","")
    new4 = new3.replace(" ","")
    numerofinal = int(new4, 2)
    return numerofinal

def RBN(K,N,T):
#    K = 2      # number of connections
#    N = 200     # number of nodes, indexed 0 .. N-1
#    T = 100    # timesteps
    Con = np.apply_along_axis(np.random.permutation, 1, np.tile(range(N), (N,1) ))[:, 0:K]
    Rules=np.random.randint(0, 2, [pow(2,K),N])
    #print(Rules)
    
    State = np.zeros((T+1,N),dtype=int)
    State[0] = np.random.randint(0, 2, N)
    for t in range(T):  # 0 .. T-1
        p=State[t,Con]
        for n in range(N):
            lista=p.tolist()
            p3=conv(lista[n])
            State[t+1][n]=Rules[p3][n]
    
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.show()
    
    atractores,indices,ContAtrac=np.unique(State,return_index=True, return_counts=True,axis=0)
    
    #print(atractores)
    #    print("Indices "+str(indices))
    #    print("Repeticiones "+ str(ContAtrac))
    
    Atr=np.size(np.where(ContAtrac>1))
    
    print("Atractores: "+str(Atr))    
    
    #print(np.argmax(np.bincount(State[:,0])))
    #    return np.argmax(np.bincount(State[:,0]))
    
    
    return np.argmax(np.bincount(State[:,0]))