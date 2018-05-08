import time
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import hamming
import multiprocessing
from functools import partial

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
             self.Bool = np.random.choice([0, 1], size=(N, 2**self.K), p=[1-p, p]) # N random boolean functions, a list of 2^k  ones and zeros.
             #self.Bool = np.random.randint(0, 2, size=(N, 2**self.K))  # N random boolean functions, a list of 2^k  ones and zeros.
        else:
            Kv=np.random.poisson(self.K, N)
            Kv[np.where(Kv>N)]=N
            #print(np.mean(Kv))
            maximo=np.amax(Kv)
            
            self.Con=np.zeros((N+1, maximo),dtype=int)
            self.Bool=np.zeros((N+1, 2**maximo),dtype=int)
            
            for i in range(N):
                self.Con[i+1, 0:Kv[i]] = np.random.choice(N, Kv[i], replace=False)+1
                self.Bool[i+1, 0:2**Kv[i]] = (np.random.choice([0, 1], size=2**Kv[i], p=[1-p, p]))
        return
    
    def RunNet(self, T, initial=[], M=0, O=0):
        """
        Con= matrix of connections
        Bool= lookup table
        T = timesteps
        initial = initial state (random if empty)
        M = how many perturbations
        O = how often the perturbations take place
        """
        Pow = 2**np.arange(np.size(self.Con, 1)) # [ 1 2 4 ... ], for converting inputs to numerical value
        
        if(type(self.K) is int):
            State = np.zeros((T+1,self.N),dtype=int)
            if np.array_equal(initial, []):
                State[0] = np.random.randint(0, 2, self.N) 
            else:
                State[0] = initial
        else:
            State = np.zeros((T+1,self.N+1),dtype=int)
            if np.array_equal(initial, []):
                State[0] = np.append([0], np.random.randint(0, 2, self.N))
            else:
                State[0] = np.append([0],initial)
            self.Bool[np.where(self.Con[:,0]==0),0] = State[0, np.where(self.Con[:,0]==0)] # if node doesn't have conections not change 
        
        for t in range(T):  # 0 .. T-1
                State[t+1] = self.Bool[:, np.sum(Pow * State[t,self.Con],1)].diagonal()
                if ( M and O ) != 0:  #Perturbations 
                    if t%O == 0:
                        State[t+1,  np.random.choice(self.N, size=M, replace=False)] = np.random.randint(0, 2, M)
                    
                
        if(type(self.K) is int):
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
                
                State=red.RunNet(T, initial)
                unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0)      
                A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion
                
                if not(A.tolist() in attractList):  #if A is not in attractList then add it
                    attractList.append(A.tolist())
                
        else:
            for i in range(runs):
                State=red.RunNet(T)
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
    
    def antifragile(self, T, runs=1):
        """
        plot antifragility of RBN
        """
        
        f=np.zeros(int(self.N/2))
        O=1
        pool = multiprocessing.Pool()
        
        for j in range(runs):
            initial = np.random.randint(0, 2, self.N)
            State=self.RunNet(T, initial)
            C0 = entropia(State)
            f+=pool.map(partial(self.func, T=T, initial=initial, O=O, C0=C0), range(1, int(self.N/2)+1))
        f/=runs # average fragility by perturbation
        plt.plot(f)
        plt.ylabel("Fragility")
        plt.xlabel("Perturbations")
        plt.show()
        
        return np.amin(f), np.argmin(f) # maximum antifragility
    
    def func(self, i, T, initial, O, C0):
        f=np.zeros(int(self.N/2))
        State=self.RunNet(T, initial, i, O)
        C = entropia(State)
        f=fragility(C, C0, i, O, self.N, T)
        return f

def entropia(state):
    """
    Measuring Complexity Based on Shanon Entropy 
    state = matrix of a RBN states
    """
    p1=np.sum(state, axis=0)/np.size(state, 0)
    p0=1-p1
    np.place(p0, p0==0, 1)
    np.place(p1, p1==0, 1)
    
    #column by column
    E=-(p0*np.log2(p0)+p1*np.log2(p1)) #Shannon Entropy
    C=4*E*(1-E) #Complexity
    return C

def fragility(C, C0, M, O, N, T):
    dx =(M*(T/O))/(N*T) # degree of perturbation
    sigma = np.mean(C-C0) # degree of satisfaction
    return -sigma*dx

def plots(red, T, M, O):
    initial = np.random.randint(0, 2, red.N)
    State=red.RunNet(T, initial)
    plt.imshow(State, cmap='Greys', interpolation='None', aspect='auto')
    plt.show()
    C1 = entropia(State)
    print("Complejidad: "+str(np.mean(C1)))
    initial2 = np.random.randint(0, 2, red.N)
    State2=red.RunNet(T, initial2)
    plt.imshow(State2, cmap='Greys', interpolation='None', aspect='auto')
    plt.show()
    C3 = entropia(State2)
    print("Complejidad: "+str(np.mean(C3)))
    print("\nDistancia inicial: ")
    print(hamming(State[0], State2[0])) #distance
    print("Distancia final: ")
    print(hamming(State[T-1], State2[T-1])) #distance
    
    State=red.RunNet(T, initial, M, O)
    plt.imshow(State, cmap='Greys', interpolation='None', aspect='auto')
    plt.show()
    C2 = entropia(State)
    print("Complejidad: "+str(np.mean(C2)))
    print("Fragilidad: "+str(fragility(C2, C1, M, O, red.N, T))) 
    State2=red.RunNet(T, initial2, M, O)
    plt.imshow(State2, cmap='Greys', interpolation='None', aspect='auto')
    plt.show()
    C4 = entropia(State2)
    print("Complejidad: "+str(np.mean(C4)))
    print("Fragilidad: "+str(fragility(C4, C3, M, O, red.N, T)))
    print("\nDistancia inicial: ")
    print(hamming(State[0], State2[0])) #distance
    print("Distancia final: ")
    print(hamming(State[T-1], State2[T-1])) #distance
    return
    
if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    
    start_time = time.time()
    
    K=2.0
    N=100
    p=0.5
    T=200
    
    M=5 # how many perturbations
    O=1 # how often the perturbations take place
    
    red=RBN()
    red.CreateNet(K, N, p)
#    print(red.Con)
#    print(red.Bool)
#    red.RBNSort()
    
    plots(red, T, M, O)
    
    print(red.antifragile(T, runs=100))
    
    A=red.Attractors(T, runs=1000)
    print("\nAttractores: ")
    print(len(A))
    edos=0
    for x in A:
        edos+=len(x)
    print("Longitud promedio de Attractores: ")
    print(edos/len(A))
    if(edos!= 0):
        print(str(len(A)/(edos)*100)+"%")
    
    print("--- %s seconds ---" % (time.time() - start_time))