import matplotlib.pyplot as plt
import numpy as np
import scipy

class RBN:    
    def CreateNet(self, K, N, p):
        """
        K = number of connections
        N = number of nodes, indexed 0 .. N-1
        """
        
        Pow = 2**np.arange(K) # [ 1 2 4 ... ], for converting inputs to numerical value
        Pow=np.flip(Pow, 0)
        self.Con = np.apply_along_axis(np.random.permutation, 1, np.tile(range(N), (N,1) ))[:, 0:K]
        self.Bool = np.random.choice([False, True], size=(N, 2**K), p=[1-p, p])
        #RBNSort(N, K)
        
        #np.savetxt('test.txt', Bool,  delimiter=',', fmt='%i', newline="\n")
        
    def RunNet(self, T, initial=[]):
        """
        Con= matrix of connections
        Bool= lookup table
        T = timesteps
        """
        
        K=self.Con[0].size
        N=len(self.Con)
        
        Pow = 2**np.arange(K) # [ 1 2 4 ... ], for converting inputs to numerical value
        Pow=np.flip(Pow, 0)
        
        State = np.zeros((T+1,N),dtype=bool)
        
        if np.array_equal(initial, []):
            State[0] = np.random.randint(0, 2, N) 
        else:
            State[0] = initial
        for t in range(T):  # 0 .. T-1
            State[t+1] = self.Bool[:, np.sum(Pow * State[t,self.Con],1)].diagonal()
        
        return(State)
        
    def AttractorsRand(self, N):
        attractList=[]
        for i in range(N):
            State=red.RunNet(T)
            unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0)      
            A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion
            
            agregar=True
            
            for l in attractList:
                if np.array_equal(l, A):
                    agregar=False
            if agregar:
                attractList.append(A)
                
        return attractList
    
    def Attractors(self):
        attractList=[]
        N=len(self.Con)
        
        for i in range(np.power(2,N)):
            
            initial=[x=='1' for x in format(i, '0'+str(N)+'b')]
            
            State=red.RunNet(T, initial)
            unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0)      
            A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion
            
            agregar=True
            
            for l in attractList:
                if np.array_equal(l, A):
                    agregar=False
            if agregar:
                attractList.append(A)
                
        return attractList
        
    def RBNSort(self, N, K):
        SRun = 5     # sorting runs
        ST = 200     # sorting timesteps
        State = np.zeros((ST+1,N),dtype=bool)
        Totals = np.zeros((N),dtype=int)
        Pow = 2**np.arange(K) # [ 1 2 4 ... ], for converting inputs to numerical value
        Pow=np.flip(Pow, 0)
        
        for r in range(SRun):
            for t in range(ST):
                State[t+1] = self.Bool[:, np.sum(Pow * State[t,self.Con],1)].diagonal()
                Totals = Totals + State[t+1]
            State[0] = np.random.randint(0, 2, N) # new initial random state
            
        Index = np.argsort(Totals)    # permutation indexes for sorted order
        self.Bool = self.Bool[Index]         # permute the boolean functions
        self.Con = self.Con[Index]           # permute the connections
        
        InvIndex = np.argsort(Index)  # inverse permutation
        self.Con = InvIndex[self.Con]        # relabel the connections
        return [self.Con, self.Bool]
    


    
if __name__ == '__main__':
    import time
    start_time = time.time()
    
    K=2
    N=500
    p=0.5
    T=200
    
    red=RBN()
    red.CreateNet(K, N, p)
    
    initial = np.zeros(N,dtype=bool)
    State=red.RunNet(T,  initial)
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.show()
    
    State2=red.RunNet(T)
    plt.imshow(State2, cmap='Greys', interpolation='None')
    plt.show()
    
    print("Distancia inicial: ")
    print(scipy.spatial.distance.hamming(State[0], State2[0]))
    print("Distancia final: ")
    print(scipy.spatial.distance.hamming(State[T-1], State2[T-1]))
    
    A=red.AttractorsRand(1000)
    print("Attractors: ")
    print(len(A))
    #print(1.*A)
    edos=0
    for x in A:
        edos+=x.size
    edos/=N
    print(str(len(A)/(edos)*100)+"%")
    
    print("--- %s seconds ---" % (time.time() - start_time))