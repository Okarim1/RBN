# -*- coding: utf-8 -*-
"""
@author: okarim
"""
import time
import numpy as np
import multiprocessing
from functools import partial
import matplotlib.pyplot as plt
import pandas as pd

class RBN:
    def CreateNet(self, K, N, p):
        """
        K = number of connections
        N = number of nodes, indexed 0 .. N-1
        p = probability of one
        """
        self.K=K
        self.N=N
        if(type(self.K) is int):
             self.Con = np.apply_along_axis(np.random.permutation, 1, np.tile(range(N), (N,1) ))[:, 0:self.K]
             #self.Bool = np.random.choice([0, 1], size=(N, 2**self.K), p=[1-p, p]) # N random Boolean functions, a list of 2^K ones and zeros.
             self.Bool = np.random.randint(0, 2, size=(N, 2**self.K))  # N random boolean functions, a list of 2^k  ones and zeros.
        else:
            Kv=np.random.poisson(self.K, N)
            Kv=Kv/(np.mean(Kv)/K)
            Kv=np.ceil(Kv)
            Kv=Kv.astype(int)
            Kv[np.where(Kv>N)]=N
            maximo=np.amax(Kv)
            
            self.Con=np.zeros((N+1, maximo),dtype=int)
            self.Bool=np.zeros((N+1, 2**maximo),dtype=int)
            
            for i in range(N):
                self.Con[i+1, 0:Kv[i]] = np.random.choice(N, Kv[i], replace=False)+1
                self.Bool[i+1, 0:2**Kv[i]] = (np.random.choice([0, 1], size=2**Kv[i], p=[1-p, p]))
        return
    
    def CreateBioNet(self, b=1):
        if(b==1):
            self.N=18
            
            data=pd.read_csv("Con_CD4+TCELL.csv", sep=",", header=-1)
            data=data.fillna(0)
            x=np.array(data.loc[:, 1:]).astype(int)
            self.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x])
            
            data=pd.read_csv("Bool_CD4+TCELL.csv", sep=",", header=0)
            data=data.fillna(0)
            y=np.array(data.loc[:, :]).astype(int).transpose()
            self.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y])
            self.e=6
            self.K=np.count_nonzero(self.Con)/(self.N-self.e)
            
        elif(b==2):
            self.N=20
            
            data=pd.read_csv("Con_Mammalian Cell Cycle.csv", sep=",", header=-1)
            data=data.fillna(0)
            x=np.array(data.loc[:, 1:]).astype(int)
            self.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x])
            
            data=pd.read_csv("Bool_Mammalian Cell Cycle.csv", sep=",", header=0)
            data=data.fillna(0)
            y=np.array(data.loc[:, :]).astype(int).transpose()
            self.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y])
            self.e=1
            self.K=np.count_nonzero(self.Con)/(self.N-self.e)
            
        elif(b==3):
            self.N=15
            
            data=pd.read_csv("Con_Cardiac development.csv", sep=",", header=-1)
            data=data.fillna(0)
            x=np.array(data.loc[:, 1:]).astype(int)
            self.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x])
            
            data=pd.read_csv("Bool_Cardiac development.csv", sep=",", header=0)
            data=data.fillna(0)
            y=np.array(data.loc[:, :]).astype(int).transpose()
            self.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y])
            self.e=1
            self.K=np.count_nonzero(self.Con)/(self.N-self.e)
            
        elif(b==4):
            self.N=12
            
            data=pd.read_csv("Con_Metabolic.csv", sep=",", header=-1)
            data=data.fillna(0)
            x=np.array(data.loc[:, 1:]).astype(int)
            self.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x])
            
            data=pd.read_csv("Bool_Metabolic.csv", sep=",", header=0)
            data=data.fillna(0)
            y=np.array(data.loc[:, :]).astype(int).transpose()
            self.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y])
            self.e=1
            self.K=np.count_nonzero(self.Con)/(self.N-self.e)
            
        elif(b==5):
            self.N=28
            
            data=pd.read_csv("Con_Death.csv", sep=",", header=-1)
            data=data.fillna(0)
            x=np.array(data.loc[:, 1:]).astype(int)
            self.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x])
            
            data=pd.read_csv("Bool_Death.csv", sep=",", header=0)
            data=data.fillna(0)
            y=np.array(data.loc[:, :]).astype(int).transpose()
            self.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y])
            self.e=3
            self.K=np.count_nonzero(self.Con)/(self.N-self.e)
            
        elif(b==6):
            self.N=14
            
            data=pd.read_csv("Con_Arabidopsis.csv", sep=",", header=-1)
            data=data.fillna(0)
            x=np.array(data.loc[:, 1:]).astype(int)
            self.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x])
            
            data=pd.read_csv("Bool_Arabidopsis.csv", sep=",", header=0)
            data=data.fillna(0)
            y=np.array(data.loc[:, :]).astype(int).transpose()
            self.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y])
            self.e=0
            self.K=np.count_nonzero(self.Con)/(self.N-self.e)
            
        elif(b==7):
            self.N=32
            
            data=pd.read_csv("Con_Tumour Cell.csv", sep=",", header=-1)
            data=data.fillna(0)
            x=np.array(data.loc[:, 1:]).astype(int)
            self.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x])
            
            data=pd.read_csv("Bool_Tumour Cell.csv", sep=",", header=0)
            data=data.fillna(0)
            y=np.array(data.loc[:, :]).astype(int).transpose()
            self.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y])
            self.e=2
            self.K=np.count_nonzero(self.Con)/(self.N-self.e)
        
        return 
    
    def RunNet(self, T, initial=[], X=0, O=0, Bio=False):
        """
        Con= matrix of connections
        Bool= lookup table
        T = timesteps
        initial = initial state (random if empty)
        X = how many perturbations
        O = how often the perturbations take place
        """
        Pow = 2**np.arange(np.size(self.Con, 1)) # [ 1 2 4 ... ], for converting inputs to numerical value
        
        if(type(self.K) is int):
            a=0
            State = np.zeros((T+1,self.N),dtype=int)
            if np.array_equal(initial, []):
                State[0] = np.random.randint(0, 2, self.N) 
            else:
                State[0] = initial
        else:
            a=1
            State = np.zeros((T+1,self.N+1),dtype=int)
            if np.array_equal(initial, []):
                State[0] = np.append([0], np.random.randint(0, 2, self.N))
            else:
                State[0] = np.append([0],initial)
            self.Bool[np.where(self.Con[:,0]==0),0] = State[0, np.where(self.Con[:,0]==0)] # if node doesn't have conections not change 
        
        for t in range(T):  # 0 .. T-1
                self.Bool[np.where(self.Con[:,0]==0),0] = State[t, np.where(self.Con[:,0]==0)]
                State[t+1] = self.Bool[:, np.sum(Pow * State[t,self.Con],1)].diagonal()
                if ( X and O ) != 0:  #Perturbations 
                    if t%O == 0:
                        State[t+1,  np.random.choice(self.N, size=X, replace=False)+a] = np.random.randint(0, 2, X)
                if(Bio and self.e>0):
                    State[t+1, self.N-self.e+1:]=np.random.randint(0, 2, self.e)
        if(Bio and self.e>0):
            return(State[:,1:-self.e])
        elif(type(self.K) is int):
            return(State)
        else:
            return(State[:,1:])
    
    def Attractors(self, T, runs=0):
        """
        List of Attractors of R random initial states
        runs = number of runs (if 0 then List of Attractors of every possible initial state)
        T = timesteps
        """
        attractList=[]
        if runs == 0 :
            for i in range(np.power(2,self.N)):
                initial=[x=='1' for x in format(i, '0'+str(self.N)+'b')]
                
                State=self.RunNet(T, initial)
                unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0)      
                A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion
                
                if not(A.tolist() in attractList):  #if A is not in attractList then add it
                    attractList.append(A.tolist())
                
        else:
            for i in range(runs):
                State=self.RunNet(T)
                unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0)      
                A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion
    
                if not(A.tolist() in attractList):  #if A is not in attractList then add it
                    attractList.append(A.tolist())
                
        return attractList
        
    def RBNSort(self): 
        """
        Sort the nodes by their overall activity
        """
        SRun = 5     # sorting runs
        ST = 200     # sorting timesteps
        Totals = np.zeros(self.N,dtype=int)
        
        for r in range(SRun):
            State=self.RunNet(ST)
            Totals = Totals + np.sum(State, 0)
        
        Index = np.argsort(Totals)    # permutation indexes for sorted order
        
        if(type(self.K) is int):
            self.Bool = self.Bool[Index]         # permute the boolean functions
            self.Con = self.Con[Index]           # permute the connections
            
            InvIndex = np.argsort(Index)         # inverse permutation
            self.Con = InvIndex[self.Con]        # relabel the connections
        else:
            self.Bool[1:] = self.Bool[Index+1]         # permute the boolean functions
            self.Con[1:] = self.Con[Index+1]           # permute the connections
            InvIndex = np.append([-1], np.argsort(Index)) # inverse permutation
            self.Con[1:] = InvIndex[self.Con[1:]]+1        # relabel the connections
        return
    
    def antifragile(self, T, runs=1, X=None, O=None, fraction=1):
        """
        plot antifragility of RBN
        """
        
        f=np.zeros(int(self.N/fraction))
        pool = multiprocessing.Pool()
        
        for j in range(runs):
            initial = np.random.randint(0, 2, self.N)
            State=self.RunNet(2*T, initial)
            C0 = complexity(State[-T:])
            if(O!=None):
                f+=pool.map(partial(self.func, T=T, initial=initial, O=O, C0=C0, fraction=fraction), range(1, int(self.N/fraction)+1))
            elif(X!=None):
                f+=pool.map(partial(self.func2, T=T, initial=initial, X=X, C0=C0, fraction=fraction), range(1, int(T/fraction)+1))
        f/=runs # average fragility by perturbation
        pool.close()
        return f
    
    def func2(self, i, T, initial, X, C0, fraction=1):
        f=np.zeros(int(self.N/fraction))
        State=self.RunNet(T*2, initial, X, i)
        C = complexity(State[-T:])
        f=fragility(C, C0, X, i, self.N, T)
        return f
    
    def func(self, X, T, initial, O, C0, fraction=1):
        f=np.zeros(int(self.N/fraction))
        State=self.RunNet(T*2, initial, X, O)
        C = complexity(State[-T:])
        f=fragility(C, C0, X, O, self.N, T)
        return f
#        if f < 0:
#            return 1
#        else:
#            return 0
        

def complexity(state):
    """
    Measuring Complexity Based on Shanon Entropy 
    state = matrix of a RBN states
    """
    p1=np.sum(state, axis=0)/np.size(state, 0)
    p0=1-p1
    np.place(p0, p0==0.0, 1.0)
    np.place(p1, p1==0.0, 1.0)
    #column by column
    E=-(p0*np.log2(p0)+p1*np.log2(p1)) #Shannon Entropy
    E=np.mean(E)
    C=4*E*(1-E) #Complexity
    return C

def fragility(C, C0, X, O, N, T):
    """
    C0 = initial complexity
    C = final complexity
    X = how many perturbations
    O = how often the perturbations take place
    N = number of nodes, indexed 0 .. N-1
    T = timesteps
    """
    dx =(X*(T/O))/(N*T) # degree of perturbation
    sigma = np.mean(C-C0) # degree of satisfaction
    #return sigma
    return -sigma*dx
    
if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    
    start_time = time.time()
    
    K=3
    N=20
    p=0.5
    T=20
    X=2
    O=1
    
    red=RBN()
    red.CreateNet(K, N, p)
#    print(red.Con)
#    print(red.Bool)
#    
#    A=red.Attractors(T, runs=1000)
#    print("\nAttractores: ")
#    print(len(A))
#    edos=0
#    for i in A:
#        edos+=len(i)
#    print("Longitud promedio de Attractores: ")
#    print(edos/len(A))
#    if(edos!= 0):
#        print(str(len(A)/(edos)*100)+"%")
    
    initial=np.random.randint(0, 2, N)
    
    State=red.RunNet(2*T, initial)
    fig= plt.figure()
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.xlabel('Node')
    plt.ylabel('Time')
    plt.title("Without perturbations")
    
    plt.savefig("Figures/Figure_1b.eps")
    plt.show()
    
    
    
    C0=complexity(State[-T:])
    
    print(C0)
    
    State=red.RunNet(2*T, initial, X=X, O=O)
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.xlabel('Node')
    plt.ylabel('Time')
    plt.title("With perturbations")
    plt.savefig("Figures/Figure_1b2.eps")
    plt.show()
    
    C=complexity(State[-T:])
    
    print(C)
    
    print(fragility(C,C0,X,O,N,T))
        
    print("--- %s seconds ---" % (time.time() - start_time))