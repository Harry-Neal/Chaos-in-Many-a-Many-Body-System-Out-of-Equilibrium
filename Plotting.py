import numpy as np
import matplotlib.pyplot as plt

f_name = "Av.dat" #enter file name
data = np.loadtxt(f_name)
x = data[:,0]
y = data[:,1]


plt.plot(x,y)
plt.xlabel('$\\tau$')
plt.ylabel('$Mean Energy Density$')
plt.show()