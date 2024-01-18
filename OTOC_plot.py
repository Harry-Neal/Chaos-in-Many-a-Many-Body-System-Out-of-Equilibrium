import matplotlib.pyplot as plt
import numpy as np
import scipy as scpy

OTOC = np.loadtxt('run/OTOC_good_equlib/OTOC_avg.dat')
params = np.loadtxt('run/OTOC_good_equlib/plotting_params_OTOC.dat')

t_max,x_max = np.shape(OTOC)
x = np.arange(0,x_max,1)
t = np.arange(0,t_max,params[1])

#fig,ax1 = plt.subplots(1)
#times = np.array([10,20,30,40,50,60,70,80,90,100])
#for time in times:
#    t_sample = int(np.argwhere(t == time))
#    ax1.plot(x[0:200],OTOC[t_sample,0:200],label=f't={t_sample}')

#ax1.set_xlabel('x')
#ax1.set_ylabel('D(x,t)')
#ax1.legend()

#fig,ax2 = plt.subplots(1)
#times = np.array([10,20,30,40,50,60,70,80,90,100])
#y=np.log(OTOC[:,1:200]/(0.05**2))

#def f(x,v,mu):
#    return mu*(1 - (v*(x))**2)


#for time in times:
#    t_sample = int(np.argwhere(t == time)) 
#    ax2.plot(x[1:200]/time,y[t_sample,:]/(2*time),label=f't={time}')

#ax2.set_xlabel('x/t')
#ax2.set_ylabel('$ln(D(x,t)/ \epsilon^2 )/2t$')
#ax2.legend()

OTOC = np.roll(OTOC,x_max//2,axis=1)
x = np.arange(-x_max//2,x_max//2,1)
fig,ax3 = plt.subplots(1)
plt.imshow(OTOC,extent=[min(x), max(x), min(t), max(t)],origin='lower',cmap='bwr')
ax3.set_aspect('equal')
ax3.set_xlabel('x')
ax3.set_ylabel('t')
ax3.plot(x,(1/1.6417)*np.abs(x)+5,'black')
cax = plt.axes((0.85, 0.1, 0.05, 0.8))
plt.colorbar(cax=cax)
cax.set_ylabel('D(x,t)')
plt.show()
