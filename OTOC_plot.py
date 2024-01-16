import matplotlib.pyplot as plt
import numpy as np
OTOC = np.loadtxt('OTOC_avg.dat')
params = np.loadtxt('plotting_params_OTOC.dat')

t_max,x_max = np.shape(OTOC)
x = np.arange(0,x_max,1)
t = np.arange(0,t_max,params[1])

fig,ax1 = plt.subplots(1)
times = np.array([10,20,30,40,50,60,70,80,90,100])
for time in times:
    t_sample = int(np.argwhere(t == time))
    ax1.plot(x[0:400],OTOC[t_sample,0:400],label=f't={t_sample}')

ax1.set_xlabel('x')
ax1.set_ylabel('D(x,t)')
ax1.legend()

fig,ax2 = plt.subplots(1)
times = np.array([10,20,30,40,50,60,70,80,90,100])
y=np.log(OTOC/(0.05**2))
for time in times:
    pass
    t_sample = int(np.argwhere(t == time)) 
    ax2.plot(x[0:600]/time,y[t_sample,0:600]/(2*time),label=f't={time}')

v = np.arange(0,1.8,0.1)
ax2.plot(v,0.494*(1 - (1.6417*(v))**2),'black')

ax2.set_xlabel('x/t')
ax2.set_ylabel('$ln(D(x,t)/ \epsilon^2 )/2t$')
ax2.legend()

OTOC = np.roll(OTOC,x_max//2,axis=1)
x = np.arange(-x_max//2,x_max//2,1)
fig,ax3 = plt.subplots(1)
ax3.pcolormesh(x,t,OTOC)
ax3.set_aspect('equal')
ax3.set_xlabel('x')
ax3.set_ylabel('t')
ax3.plot(x,(1/1.6417)*np.abs(x),'black')

plt.show()
