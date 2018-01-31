import matplotlib.pyplot as plt
import numpy as np
import scipy
    
def CreateNet(K, N, p):
    """
    K = number of connections
    N = number of nodes, indexed 0 .. N-1
    """
    
    Pow = 2**np.arange(K) # [ 1 2 4 ... ], for converting inputs to numerical value
    Pow=np.flip(Pow, 0)
    Con = np.apply_along_axis(np.random.permutation, 1, np.tile(range(N), (N,1) ))[:, 0:K]
    Bool = np.random.choice([False, True], size=(N, 2**K), p=[1-p, p])
    #[Con, Bool]=RBNSort(Con, Bool, N, K)
    
    #np.savetxt('test.txt', Bool,  delimiter=',', fmt='%i', newline="\n")
    return [Con, Bool]
    
def RunNet(Con, Bool, T):
    """
    Con= matrix of connections
    Bool= lookup table
    T = timesteps
    """
    
    K=Con[0].size
    N=len(Con)
    
    Pow = 2**np.arange(K) # [ 1 2 4 ... ], for converting inputs to numerical value
    Pow=np.flip(Pow, 0)
    
    State = np.zeros((T+1,N),dtype=bool)
    State[0] = np.random.randint(0, 2, N)
    for t in range(T):  # 0 .. T-1
        State[t+1] = Bool[:, np.sum(Pow * State[t,Con],1)].diagonal()
    
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.show()
    
    return(State)
    
def Attractors(State): 
    unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0)      
    A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion
            
    return A
            
    
def RBNSort(Con, Bool, N, K):
    SRun = 5     # sorting runs
    ST = 200     # sorting timesteps
    State = np.zeros((ST+1,N),dtype=bool)
    Totals = np.zeros((N),dtype=int)
    Pow = 2**np.arange(K) # [ 1 2 4 ... ], for converting inputs to numerical value
    Pow=np.flip(Pow, 0)
    
    for r in range(SRun):
        for t in range(ST):
            State[t+1] = Bool[:, np.sum(Pow * State[t,Con],1)].diagonal()
            Totals = Totals + State[t+1]
        State[0] = np.random.randint(0, 2, N) # new initial random state
        
    Index = np.argsort(Totals)    # permutation indexes for sorted order
    Bool = Bool[Index]         # permute the boolean functions
    Con = Con[Index]           # permute the connections
    
    InvIndex = np.argsort(Index)  # inverse permutation
    Con = InvIndex[Con]        # relabel the connections
    return [Con, Bool]
    
if __name__ == '__main__':
    #import timeit
    #print(timeit.timeit("RBN(2, 5, 10)", setup="from __main__ import RBN"))
    
    K=2
    N=500
    p=0.5
    T=200
    
    [Con,Bool]=CreateNet(K, N, p)
    
    State=RunNet(Con, Bool, T)
    State2=RunNet(Con, Bool, T)
    
    print("Distancia inicial: ")
    print(scipy.spatial.distance.hamming(State[0], State2[0]))
    print("Distancia final: ")
    print(scipy.spatial.distance.hamming(State[T-1], State2[T-1]))
    
    
    A=Attractors(State)
    print("Attractors: ")
    print(1.*A)
    print(len(A))
    print(str(1/len(A)*100)+"%")
    