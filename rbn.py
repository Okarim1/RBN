import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import hamming

def bool2int(x):
    y = 0
    for i,j in enumerate(x):
        y += j<<i
    return y

class RBN:    
    def CreateNet(self, K, N, p):
        """
        K = number of connections
        N = number of nodes, indexed 0 .. N-1
        """
        if(type(K) is int):
             self.Con = np.apply_along_axis(np.random.permutation, 1, np.tile(range(N), (N,1) ))[:, 0:K]
             self.Bool = np.random.choice([0, 1], size=(N, 2**K), p=[1-p, p]) # N random boolean functions, a list of 2^k  ones and zeros.
             #self.Bool = np.random.randint(0, 2, size=(N, 2**K))  # N random boolean functions, a list of 2^k  ones and zeros.
        else:
            Kv=np.random.poisson(K, N)
            #print(np.mean(Kv))
            #maximo=np.amax(Kv)
            
            #self.Con=np.zeros((maximo, N),dtype=int)
            #self.Bool=np.zeros((2**maximo, N),dtype=int)
            
            self.Con=[]
            self.Bool=[]
            for i in range(N):
                self.Con.append(np.random.choice(N, Kv[i], replace=False))
                self.Bool.append(np.random.choice([0, 1], size=2**Kv[i], p=[1-p, p]))
        
        
    def RunNet(self, T, initial=[]):
        """
        Con= matrix of connections
        Bool= lookup table
        T = timesteps
        initial = initial state (random if empty)
        """
        
        N=len(self.Con)
        
        State = np.zeros((T+1,N),dtype=int)
        
        Pow = 2**np.arange(K) # [ 1 2 4 ... ], for converting inputs to numerical value
        Pow=np.flip(Pow, 0)
        
        if np.array_equal(initial, []):
            State[0] = np.random.randint(0, 2, N) 
        else:
            State[0] = initial
        for t in range(T):  # 0 .. T-1
            if(type(self.Bool) is list):
                State[t+1]=[self.Bool[i][bool2int(State[t,self.Con[i]])] for i in range(N)]
            else:
                State[t+1] = self.Bool[:, np.sum(Pow * State[t,self.Con],1)].diagonal()
        return(State)
        
    def AttractorsRand(self, N):
        attractList=[]
        for i in range(N):
            State=red.RunNet(T)
            unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0)      
            A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion

            if not(A.tolist() in attractList):
                attractList.append(A.tolist())
                
        return attractList
    
    def Attractors(self):
        attractList=[]
        N=len(self.Con)
        
        for i in range(np.power(2,N)):
            initial=[x=='1' for x in format(i, '0'+str(N)+'b')]
            
            State=red.RunNet(T, initial)
            unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0)      
            A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion
            
            if not(A.tolist() in attractList):
                attractList.append(A.tolist())
                
        return attractList
        
    def RBNSort(self, N, K):
        SRun = 5     # sorting runs
        ST = 200     # sorting timesteps
        State = np.zeros((ST+1,N),dtype=int)
        Totals = State[0]
        #        Totals = np.zeros((N),dtype=int)
        Pow = 2**np.arange(K) # [ 1 2 4 8 16... ], for converting inputs to numerical value
        Pow=np.flip(Pow, 0)
        
        for r in range(SRun):
            for t in range(ST):
                if(type(self.Bool) is list):
                    State[t+1]=[self.Bool[i][bool2int(State[t,self.Con[i]])] for i in range(N)]
                else:
                    State[t+1] = self.Bool[:, np.sum(Pow * State[t,self.Con],1)].diagonal()
                Totals = Totals + State[t+1]
            State[0] = np.random.randint(0, 2, N) # new initial random state
   
        Index = np.argsort(Totals)    # permutation indexes for sorted order
        self.Bool = [self.Bool[i] for i in Index]   # permute the boolean functions
        self.Con = [self.Con[i] for i in Index]     # permute the connections
        
        InvIndex = np.argsort(Index)  # inverse permutation
        self.Con = [InvIndex[i] for i in self.Con]        # relabel the connections
        return [self.Con, self.Bool]
    

def entropia(state):
    p1=np.sum(state, axis=0)/np.size(state, 0)
    p0=1-p1
    np.place(p0, p0==0, 1)
    np.place(p1, p1==0, 1)
    
    E=-(p0*np.log2(p0)+p1*np.log2(p1)) #Entropia
    print("Entropia: "+str(np.mean(E)))
    
    C=4*E*(1-E) #Complejidad
    
    print("Complejidad: "+str(np.mean(C)))

    
if __name__ == '__main__':
    import time
    start_time = time.time()
    
    K=2.5
    N=100
    p=0.5
    T=500
    
    red=RBN()
    red.CreateNet(K, N, p)
#    print(red.Con)
#    print(red.Bool)
#    red.RBNSort(N, K)
    
    initial = np.zeros(N,dtype=int)
    State=red.RunNet(T,  initial)
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.show()
    
    entropia(State)
    
    State2=red.RunNet(T)
    plt.imshow(State2, cmap='Greys', interpolation='None')
    plt.show()
    
    entropia(State2)
    
    print("Distancia inicial: ")
    print(hamming(State[0], State2[0])) #distancia
    print("Distancia final: ")
    print(hamming(State[T-1], State2[T-1])) #distancia
    
    A=red.AttractorsRand(100)
    print("Attractores: ")
    print(len(A))
    
    edos=0
    for x in A:
        edos+=len(x)
    
    print("Longitud promedio de Attractores: ")
    print(edos/len(A))
    
    print(str(len(A)/(edos)*100)+"%")
    
    print("--- %s seconds ---" % (time.time() - start_time))