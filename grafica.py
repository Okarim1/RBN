import matplotlib.pyplot as plt
import numpy as np

data = np.load("valores5.npz")
rango=np.arange(0.1, 6.0, 0.1)

vfm= data['arr_0']
vpm= data['arr_1']
vfbp= data['arr_2']
vpbp= data['arr_3']

plt.plot(rango, vfm, 'o')
c=np.polyfit(np.log(rango), vfm, 1)
y=c[0]*np.log(rango)+c[1]
plt.plot(rango, y)
plt.ylabel("Fragility")
plt.xlabel("K")

plt.show()   
plt.plot(rango, vpm, 'o')
c=np.polyfit(rango, np.log(vpm), 1, w=np.sqrt(vpm))
y=np.exp(c[1])*np.exp(rango*c[0])
plt.plot(rango, y)
plt.ylabel("Perturbations")
plt.xlabel("K")

plt.show()   
plt.plot(rango, vfbp, 'o')
c=np.polyfit(np.log(rango), vfm, 1)
y=c[0]*np.log(rango)+c[1]
plt.plot(rango, y, 'r', linewidth=5.0, label=str(c[0])+' log(X) '+str(c[1]))
plt.legend()
plt.ylabel("Fragility")
plt.xlabel("K")


plt.show()   
plt.plot(rango, vfbp, 'o')
c=np.polyfit(rango, vfm, 3)
y=np.polyval(c, rango)
plt.plot(rango, y, 'r', linewidth=5.0, 
         label=str("%.4f"%c[0])+' X^3+'+str("%.4f"%c[1])+' X^2+'+str("%.4f"%c[2])+' X+'+str("%.4f"%c[3]))
plt.legend()
plt.ylabel("Fragility")
plt.xlabel("K")

plt.show()
plt.plot(rango, vpbp, 'o')
c=np.polyfit(rango, np.log(vpm), 1, w=np.sqrt(vpm))
y=np.exp(c[1])*np.exp(rango*c[0])
plt.plot(rango, y, 'r', linewidth=5.0, label=str("%.4f"%np.exp(c[1]))+' e^('+str("%.4f"%c[0])+' X)')
plt.legend()
plt.ylabel("Perturbations")
plt.xlabel("K")