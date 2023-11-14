import numpy as np
import matplotlib.pyplot as plt

f_name = "data/Av.dat" #enter file name
data = np.loadtxt(f_name)
x = data[:,0]
y = data[:,1]


plt.bar(x,y)
plt.show()