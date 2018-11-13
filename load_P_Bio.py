import matplotlib.pyplot as plt
import numpy as np
import rbn

plt.style.use('classic')

red=rbn.RBN()

red.CreateBioNet(1)
data=np.load("BioProb.npz")
Z=data['arr_0']
fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,red.N,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability for CD4+ T Cell")

red.CreateBioNet(2)
data=np.load("BioProb2.npz")
Z=data['arr_0']
fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,red.N,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability for Mammalian Cell Cycle")

red.CreateBioNet(3)
data=np.load("BioProb3.npz")
Z=data['arr_0']
fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,red.N,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability for Cardiac development")

red.CreateBioNet(4)
data=np.load("BioProb4.npz")
Z=data['arr_0']
fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,red.N,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability for Metabolic Interactions in the Gut Microbiome")

red.CreateBioNet(5)
data=np.load("BioProb5.npz")
Z=data['arr_0']
fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,red.N,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability for Death Receptor Signaling")

red.CreateBioNet(6)
data=np.load("BioProb6.npz")
Z=data['arr_0']
fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,red.N,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability for Arabidopsis thaliana Cell Cycle")

red.CreateBioNet(7)
data=np.load("BioProb7.npz")
Z=data['arr_0']
fig, ax = plt.subplots()
im = ax.imshow(Z, extent=[1,red.N,20,1])
fig.colorbar(im)
plt.xlabel('X')
plt.ylabel('O')
im.set_clim(0, 1)
plt.title("Probability for Tumour Cell Invasion and Migration")