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

data=np.load("B1.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,18,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(-0.2, 0.3)
plt.title("Fragility at K=1")

data=np.load("B2.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,18,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(-0.2, 0.3)
plt.title("Fragility at K=2")

data=np.load("B3.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,18,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(-0.2, 0.3)
plt.title("Fragility at K=3")

data=np.load("B4.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,18,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(-0.2, 0.3)
plt.title("Fragility at K=4")