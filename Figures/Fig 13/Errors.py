import matplotlib.pyplot as plt
import numpy as np

taus=np.arange(0.1,2.1,0.1)
runs=200
dts = [0.005,0.01,0.02]
markers=['x','.','d']
Av_s = np.zeros((len(taus),len(dts)))
errs = np.zeros((len(taus),len(dts)))
fig,ax=plt.subplots(1)
for j,dt in enumerate(dts):
    
    for i,tau in enumerate(taus):
        Av = np.loadtxt(f'Energy_absorbtion_error/dt={dt}/{tau:.1f}/Avg_energy.dat')
        Av_s[i,j] = Av[1]
        E  = np.loadtxt(f'Energy_absorbtion_error/dt={dt}/{tau:.1f}/energy.dat')
        errs[i,j] = (np.std(E[1]))/np.sqrt(runs)
    print(dt)
    ax.errorbar(taus,Av_s[:,j],errs[:,j],marker=markers[j],linestyle = '--',capsize=2,label=f'$\delta t = {dt:.3f}$')
ax.set_xlabel('$\\tau$',size=25)
ax.set_ylabel('e',size=25)
ax.tick_params(labelsize=18)
ax.legend(fontsize=20)
fig.set_size_inches(11,6)
plt.show()