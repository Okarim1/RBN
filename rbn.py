import matplotlib.pyplot as plt
import numpy as np
    
def RBN(K, N, T):
    """
    K = number of connections
    N = number of nodes, indexed 0 .. N-1
    T = timesteps
    """
    
    Pow = 2**np.arange(K) # [ 1 2 4 ... ], for converting inputs to numerical value
    Pow=np.flip(Pow, 0)
    Con = np.apply_along_axis(np.random.permutation, 1, np.tile(range(N), (N,1) ))[:, 0:K]
    Bool = np.random.randint(2, size=(N, 2**K), dtype=bool)
    
    #[Bool, Con]=RBNSort(Bool, Con, N, K)
    
    State = np.zeros((T+1,N),dtype=bool)
    State[0] = np.random.randint(0, 2, N)
    for t in range(T):  # 0 .. T-1
        State[t+1] = Bool[:, np.sum(Pow * State[t,Con],1)].diagonal()
    
    return(State)
    
def Attractors(State):
     
    unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0)      
    
    A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion
            
    return A
            
    
def RBNSort(Bool, Con, N, K):
    SRun = 5     # sorting runs
    ST = 200     # sorting timesteps
    State = np.zeros((ST+1,N),dtype=bool)
    Totals = State[0]
    Pow = 2**np.arange(K) # [ 1 2 4 ... ], for converting inputs to numerical value
    np.flip(Pow, 0)
    
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
    return [Bool, Con]
    
if __name__ == '__main__':
    #import timeit
    #print(timeit.timeit("RBN(2, 5, 10)", setup="from __main__ import RBN"))
    State=RBN(3, 5, 10)
        
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.show()
    
    A=Attractors(State)
    print(A)
    
    
    
    
    
    
    
    