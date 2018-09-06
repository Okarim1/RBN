import rbn
import time
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    
    start_time = time.time()
    
    K=1
    N=100
    p=0.5
    T=100
    
    M=10
    Prob=0.9
    
    X=20 # how many perturbations
    O=1 # how often the perturbations take place
    
    red=rbn.RBN()
    red.CreateNet(K, N, p)
    initial = np.random.randint(0, 2, red.N)
    State=red.RunNet(T, initial)
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.show()
    
    State=red.RunNet(T, initial, X=40, O=1)
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.show()
    
    State=red.RunNet(T, initial, X=40, O=2)
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.show()
    
    State=red.RunNet(T, initial, X=40, O=3)
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.show()
    
    State=red.RunNet(T, initial, X=40, O=4)
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.show()
    
    State=red.RunNet(T, initial, X=40, O=10)
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.show()
    
    print("--- %s seconds ---" % (time.time() - start_time))