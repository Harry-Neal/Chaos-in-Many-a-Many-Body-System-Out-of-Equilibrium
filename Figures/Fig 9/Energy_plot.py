import numpy as np
import matplotlib.pyplot as plt
path ='e_trajectories/'

fig,ax = plt.subplots(1)
colors = ['green','C0']
linestyles = ['solid','--']
i=0
for tau in [4,0.5]:
    data = np.loadtxt(path+f"tau={tau}/Av.dat")
    ax.plot(data[:,0],data[:,1],label=f'$\\tau ={tau:.1f}$',color = colors[i],linestyle=linestyles[i])
    i+=1
ax.set_xlabel('$t/\\tau$',size=25)
ax.set_ylabel('e',size=25)
ax.legend(fontsize=20)
ax.tick_params(labelsize=18)
fig.set_size_inches(11,6)
fig.tight_layout()
plt.show()