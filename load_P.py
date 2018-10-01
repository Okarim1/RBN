import matplotlib.pyplot as plt
import numpy as np

plt.style.use('classic')


data=np.load("P1.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1])
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability at K=1")

data=np.load("P2.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1])
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability at K=2")

data=np.load("P3.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1])
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability at K=3")

data=np.load("P4.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1])
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability at K=4")

data=np.load("P5.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1])
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability at K=5")