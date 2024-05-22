import numpy as np
import matplotlib.pyplot as plt

fig,ax = plt.subplots(1)
path = 'Next Nearest Neighbour/energy_absorbtion/'
djs=[0.05,0.1,0.25,0.5]
markers=['x','+','.','d']
for i,dj in enumerate(djs):
    data = np.loadtxt(path+'dj='+str(dj)+'/results.txt')
    ax.plot(data[:,0],data[:,1],label=f'$\\delta J = {dj}$',marker=markers[i])

fig.set_size_inches(11,6)
ax.set_xlabel('$\\tau$',size=25)
ax.set_ylabel('e',size=25)
ax.tick_params(labelsize=18)
ax.legend(fontsize=20)
fig.tight_layout()
plt.show()