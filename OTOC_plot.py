import matplotlib.pyplot as plt
import numpy as np
OTOC = np.loadtxt("OTOC_avg.dat")
params = np.loadtxt("plotting_params_OTOC.dat")

fig,ax = plt.subplots(1)
t_max,x_max = np.shape(OTOC)
x = np.arange(-x_max//2,x_max//2,1)
t = np.arange(0,t_max,params[1])

OTOC = np.roll(OTOC,x_max//2,axis=1)
ax.pcolormesh(x,t,OTOC)
ax.set_aspect('equal')
plt.show()