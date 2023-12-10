import numpy as np
import matplotlib.pyplot as plt

file_1 = np.loadtxt("run/200 runs dj 0.05/results.txt")
file_2 = np.loadtxt("run/200 runs dj 0.1/results.txt") 
file_3 = np.loadtxt("run/200 runs dj 0.25/results.txt")
file_4 = np.loadtxt("run/200 runs dj 0.5/results.txt")

plt.plot(file_1[:,0],file_1[:,1],color='blue',label='dj = 0.05',marker='x')
plt.plot(file_2[:,0],file_2[:,1],color='green',label='dj = 0.1',marker='x')
plt.plot(file_3[:,0],file_3[:,1],color='red',label='dj = 0.25',marker='x')
plt.plot(file_4[:,0],file_4[:,1],color='yellow',label='dj = 0.5',marker='x')

plt.xlabel('$\\tau$')
plt.ylabel('$e$')
plt.legend()
plt.show()
