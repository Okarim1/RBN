import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
#plt.style.use('classic')


data=np.load("npz/P1.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1], cmap=cm.jet)
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability at K=1")
plt.savefig("Figure_6a.eps")

data=np.load("npz/P2.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1], cmap=cm.jet)
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability at K=2")
plt.savefig("Figure_6b.eps")

data=np.load("npz/P3.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1], cmap=cm.jet)
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability at K=3")
plt.savefig("Figure_6c.eps")

data=np.load("npz/P4.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1], cmap=cm.jet)
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability at K=4")
plt.savefig("Figure_6d.eps")

data=np.load("npz/P5.npz")
Z=data['arr_0']

fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,100,30,1], cmap=cm.jet)
fig.colorbar(im, orientation='horizontal')
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability at K=5")