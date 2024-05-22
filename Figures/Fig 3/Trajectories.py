import numpy as np
import matplotlib.pyplot as plt

path1 = 'high_freq/'
Av_data1 = np.loadtxt(path1+"Av.dat")

path2 = 'low_freq/'
Av_data2 = np.loadtxt(path2+"Av.dat")

fig,ax = plt.subplots(1)

ax.plot(Av_data1[:,0],Av_data1[:,1],linestyle='--',label=f'$\\tau = 0.5$')

ax.plot(Av_data2[:,0],Av_data2[:,1],color='green',label=f'$\\tau = 4.0$')
ax.set_xlabel('Periods elapsed: $t/\\tau$',size=25)
ax.set_ylabel('e',size=25)
ax.legend(fontsize=20)
fig.set_size_inches(11,6)
fig.tight_layout()
ax.tick_params(labelsize=17)
plt.show()
