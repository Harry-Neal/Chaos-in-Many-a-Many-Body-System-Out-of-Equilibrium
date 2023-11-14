import numpy as np
import matplotlib.pyplot as plt

Analytic_dat = "data/Hist_energy_Analytic.dat" #enter file name
data1 = np.loadtxt(Analytic_dat)
x1 = data1[:,0]
y1 = data1[:,1]

Actual_dat = "data/Hist.dat"
data2 = np.loadtxt(Actual_dat)
x2 = data2[:,0]
y2 = data2[:,1]

plt.plot(x1,y1,color='r')
plt.bar(x2,y2,width=0.01,color='b')
plt.show()