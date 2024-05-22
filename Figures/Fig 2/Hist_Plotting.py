import numpy as np
import matplotlib.pyplot as plt

Analytic_dat = "Hist_Energy_Analytic.dat" #enter file name
data1 = np.loadtxt(Analytic_dat)
x1 = data1[:,0]
y1 = data1[:,1]

Single_traj = "Hist_single_trajectory.dat"
data2 = np.loadtxt(Single_traj)
x2 = data2[:,0]
y2 = data2[:,1]

ensemble = "Spin_ensemble.dat"
data3 = np.loadtxt(ensemble)

fig = plt.figure(figsize=(11,6))
plt.plot(x1,y1,color='red',label = "$\mathcal{P}_{l}(e)$")
plt.plot(x2,y2,color='green',marker='x',linestyle='--',label = "$\mathcal{P}_{T}(e)$")
plt.hist(data3,bins=np.arange(np.min(data3),np.max(data3),0.01),histtype='step',density=True,color='blue',label = "$\mathcal{P}_{\\text{can}}(e)$")
plt.xlabel('e',fontsize=25)
plt.ylabel('$\mathcal{P}(e)$',fontsize=25)
plt.xlim(-0.99,-0.33)
plt.legend(fontsize=22)
plt.tick_params(labelsize=18)
plt.show()