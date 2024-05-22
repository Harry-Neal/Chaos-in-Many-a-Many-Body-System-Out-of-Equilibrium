import numpy as np
import matplotlib.pyplot as plt
djs = [0.05,0.1,0.25,0.5]

fig = plt.figure(figsize=(11,6))
ax = plt.axes()
linestyles=['-','--','-.',':']
markers=['.','x','+','d']
for i,dj in enumerate(djs):
    file = np.loadtxt(f"Energy_absorbtion/{dj}/results.txt")
    ax.plot(file[:,0],file[:,1],label=f'$\delta j ={djs[i]}$',marker=markers[i],linestyle='--')
plt.xlabel('$\\tau$',size=25)
plt.ylabel('$e $',size=25)
plt.legend(fontsize=20)
plt.tick_params(labelsize=17)
plt.show()
