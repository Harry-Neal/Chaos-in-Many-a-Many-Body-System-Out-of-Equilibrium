import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

def fermi(x,a,b):
    return 1/(np.exp((x-a)/b) +1)

fig,ax=plt.subplots(1)

NNN=[1]
colors=['r','g','b']
for j,i in enumerate(NNN):
    folder = "Next Nearest Neighbour/Driven_NNN="+str(i)
    data = np.loadtxt(folder+"/results.txt")
    ax.plot(data[:,0],data[:,1],label=f'K = 1/{i} J',marker='x',color=colors[j])
    # ppot,pcov = opt.curve_fit(fermi,data[:,0],data[:,1])
    # ax.plot(data[:,0],fermi(data[:,0],*ppot),color=colors[j])

ax.legend()
plt.show()