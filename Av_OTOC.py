import numpy as np
import matplotlib.pyplot as plt
runs = 100

spin = np.loadtxt('OTOC_Dr_tau=0.5/Spin.dat')

time = np.arange(np.min(spin[:,0]),np.max(spin[:,0]),1)
E = np.zeros(len(time))
for i in range(runs):
    E = E + spin[0+i*len(time):len(time)+i*len(time),1]

E=E/runs
fig,ax = plt.subplots(1)
ax.plot(time,E)
ax.set_xlabel('$t/\\tau$')
ax.set_ylabel('$e$')
plt.show()