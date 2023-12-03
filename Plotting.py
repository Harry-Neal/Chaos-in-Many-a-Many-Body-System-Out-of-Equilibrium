import numpy as np
import matplotlib.pyplot as plt

f_name = "results.txt" #enter file name
data = np.loadtxt(f_name)
x = data[:,0]
y = data[:,1]


plt.plot(x,y)
plt.xlabel('$\\tau$')
plt.ylabel('$e$')
plt.show()

