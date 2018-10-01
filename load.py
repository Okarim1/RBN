import matplotlib.pyplot as plt
import numpy as np

plt.style.use('classic')

data=np.load("Bio.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,18,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(-0.2, 0.3)
plt.title("Fragility for biological Boolean network")

data=np.load("Z1.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1])
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(-0.2, 0.3)
plt.title("Fragility at K=1")

data=np.load("Z2.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1])
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(-0.2, 0.3)
plt.title("Fragility at K=2")

data=np.load("Z3.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1])
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(-0.2, 0.3)
plt.title("Fragility at K=3")

data=np.load("Z4.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1])
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(-0.2, 0.3)
plt.title("Fragility at K=4")

data=np.load("Z5.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1])
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(-0.2, 0.3)
plt.title("Fragility at K=5")