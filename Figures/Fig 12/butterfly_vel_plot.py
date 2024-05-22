import matplotlib.pyplot as plt
import numpy as np
import scipy as scpy
import scipy.optimize as opt
import re

def line(x,m,c):
    return m*x + c

djs=[0.05,0.1,0.25,0.5]
taus = np.arange(0.1,1.6,0.1)

#read OTOC data
v_bs = np.zeros((len(taus),len(djs)))
for m,dj in enumerate(djs):
    for l,tau in enumerate(taus):
        path = 'OTOC_dr/'+f'{dj}'+'/'+f'{tau:0.1f}'+'/'
        OTOC = np.loadtxt(path+'OTOC_Dr.dat')

        #read parameters file
        with open(path+'Parameters.dat') as f:
            params = f.read()

        #extract time step (tau) between points
        parameters = (re.findall(r'\w+', params))
        for i in range(len(parameters)):
            if parameters[i] == 'tau':
                tau = eval(parameters[i+1]+'.'+parameters[i+2])
            if parameters[i] =='epsilon':
                epsilon = eval(parameters[i+1]+'.'+parameters[i+2])
            if parameters[i] == 'T':
                T = eval(parameters[i+1]+'.'+parameters[i+2])
            if parameters[i] == 'dt':
                dt = eval(parameters[i+1]+'.'+parameters[i+2])

        #declare range for OTOC data
        t_max,x_max = np.shape(OTOC)
        x = np.arange(0,x_max,1)
        t = np.arange(0,T+tau,tau)
        tds = np.zeros(len(x[0:x_max//2-20]))
        threshold = 0.1
        for j,pos in enumerate(x[0:x_max//2-20]):
            for k,time in enumerate(t):
                if OTOC[k,j] >= threshold:
                        t_d = time
                        break
                else:
                        pass
            tds[j] = t_d

        ppot,pcov = opt.curve_fit(line,x[0:x_max//2-20],tds)
        sd = np.sqrt(np.diag(pcov))[0]
        v_b = 1/ppot[0]
        v_bs[l,m] = v_b

fig = plt.figure()
ax = plt.axes()
markers = ['x','+','.','d']
for i,dj in enumerate(djs):
    ax.plot(taus,v_bs[:,i],label=f'$\delta J ={dj}$',marker=markers[i],linestyle='--')

ax.set_xlabel('$\\tau$',size=25)
ax.set_ylabel('$v_b$',size=25)
ax.legend(fontsize=20)
ax.tick_params(labelsize=18)
fig.set_size_inches(11,6)
fig.tight_layout
plt.show()
