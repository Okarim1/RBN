import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import hamming

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
        
        
    def RunNet(self, T, initial=[]):
        """
        Con= matrix of connections
        Bool= lookup table
        T = timesteps
        initial = initial state (random if empty)
        """
        Pow = 2**np.arange(np.size(self.Con, 1)) # [ 1 2 4 ... ], for converting inputs to numerical value
        
        if(type(self.K) is int):
            State = np.zeros((T+1,self.N),dtype=int)
            if np.array_equal(initial, []):
                State[0] = np.random.randint(0, 2, self.N) 
            else:
                State[0] = initial
            for t in range(T):  # 0 .. T-1
                State[t+1] = self.Bool[:, np.sum(Pow * State[t,self.Con],1)].diagonal()
        else:
            State = np.zeros((T+1,self.N+1),dtype=int)
            if np.array_equal(initial, []):
                State[0] = np.append([0], np.random.randint(0, 2, self.N))
            else:
                State[0] = np.append([0],initial)
            for t in range(T):  # 0 .. T-1
                State[t+1] = self.Bool[:, np.sum(Pow * State[t,self.Con],1)].diagonal()

            State=State[:,1:]
        return(State)
        
    def AttractorsRand(self, R, T):
        """
        List of Attractors of R random initial states
        R = number of runs
        T = timesteps
        """
        attractList=[]
        for i in range(R):
            State=red.RunNet(T)
            unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0)      
            A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion

            if not(A.tolist() in attractList):  #if A is not in attractList then add it
                attractList.append(A.tolist())
                
        return attractList
    
    def Attractors(self, T):
        """
        List of Attractors of every possible initial state
        T = timesteps
        """
        attractList=[]
        
        for i in range(np.power(2,self.N)):
            initial=[x=='1' for x in format(i, '0'+str(self.N)+'b')]
            
            State=red.RunNet(T, initial)
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
    return E, C

    
if __name__ == '__main__':
    import time
    start_time = time.time()
    
    K=2.0
    N=100
    p=0.5
    T=200
    
    red=RBN()
    red.CreateNet(K, N, p)
#    print(red.Con)
#    print(red.Bool)
#    red.RBNSort()
    
    initial = np.zeros(N,dtype=int)
    State=red.RunNet(T, initial)
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.show()
    
    E, C = entropia(State)
    print("Entropia: "+str(np.mean(E)))
    print("Complejidad: "+str(np.mean(C)))
    
    State2=red.RunNet(T)
    plt.imshow(State2, cmap='Greys', interpolation='None')
    plt.show()
    
    E, C = entropia(State2)
    print("Entropia: "+str(np.mean(E)))
    print("Complejidad: "+str(np.mean(C)))
    
    print("Distancia inicial: ")
    print(hamming(State[0], State2[0])) #distance
    print("Distancia final: ")
    print(hamming(State[T-1], State2[T-1])) #distance
    
    A=red.AttractorsRand(1000, T)
    print("Attractores: ")
    print(len(A))
    
    edos=0
    for x in A:
        edos+=len(x)
    
    print("Longitud promedio de Attractores: ")
    print(edos/len(A))
    
    print(str(len(A)/(edos)*100)+"%")
    
    print("--- %s seconds ---" % (time.time() - start_time))