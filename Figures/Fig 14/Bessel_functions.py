import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
from cycler import cycler

fig,ax = plt.subplots(1)
cyc= (cycler(color=['blue','orange','green','red','purple','brown']) + cycler(linestyle=['-', '--','-', '--','-','--']))
ax.set_prop_cycle(cyc)
js=np.arange(1,7,1)
x=np.arange(0,6,0.1)
for j in js:
    I = sci.special.iv(j,x)
    ax.plot(x,I,label=f'j={j}')
ax.set_ylim(0,1)
ax.set_xlim(0,5)
ax.set_ylabel('$I_j[Jt]$',size=25)
ax.set_xlabel('Jt',size=25)
ax.tick_params(labelsize=18)
ax.legend(fontsize=20)
fig.set_size_inches(11,6)
plt.show()