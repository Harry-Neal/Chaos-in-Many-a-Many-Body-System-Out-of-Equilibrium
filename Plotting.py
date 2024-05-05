import numpy as np
import matplotlib.pyplot as plt
djs = [0.05,0.1,0.25,0.5]

fig = plt.figure(figsize=(10,8))
ax = plt.axes()
for i,dj in enumerate(djs):
    file = np.loadtxt("run/dj"+str(dj)+"/results.txt")
    ax.plot(file[:,0],file[:,1],label=f'dj = {djs[i]}',marker='x')
plt.xlabel('$\\tau$',size=30)
plt.ylabel('$e$',size=30)
plt.legend(fontsize=20)
plt.tick_params(labelsize=20)
plt.show()
