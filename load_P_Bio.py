import matplotlib.pyplot as plt
import numpy as np

plt.style.use('classic')


data=np.load("BP1.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,18,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability at K=1")

data=np.load("BP2.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,18,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability at K=2")

data=np.load("BP3.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,18,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability at K=3")

data=np.load("BP4.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,18,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability at K=4")