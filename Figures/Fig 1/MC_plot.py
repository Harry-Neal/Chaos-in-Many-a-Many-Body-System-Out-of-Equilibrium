import numpy as np
import matplotlib.pyplot as plt

#Import data
data = np.loadtxt("Energy.dat")
Av_data = np.loadtxt("Av.dat")

#set number of runs
runs = 200

#find maximum MC step
MC_max = int(np.max(Av_data[:,0]))

#plot data for individual trajectories
fig,ax = plt.subplots(1)
for i in range(runs):
    ax.plot(data[i*(MC_max+1):MC_max+i*MC_max+1,0],data[i*(MC_max+1):MC_max+i*MC_max+1,1],color='lightgrey')

#plot average data
ax.plot(Av_data[:,0],Av_data[:,1],color='black')
ax.set_xlabel('Periods elapsed: $t/\\tau$',size=25)
ax.set_ylabel('e',size=25)
fig.set_size_inches(11,6)
fig.tight_layout()
plt.tick_params(labelsize=18)
plt.show()
