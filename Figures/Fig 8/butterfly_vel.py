import numpy as np
import re
import matplotlib.pyplot as plt
import scipy.optimize as opt

def line(x,m,c):
    return m*x + c

taus = np.arange(0.1,2.1,0.1)
djs=[0.05,0.1,0.25,0.35,0.5]
v_bs=np.zeros((len(taus),len(djs)))
fig = plt.figure(figsize=(10,8))
ax = plt.axes()
markers=['x','+','.','d','s']
for m,dj in enumerate(djs):
    for l,tau_cur in enumerate(taus):
        path = 'OTOC_driving/'+ str(dj) +'/' + f'{tau_cur:.1f}' + '/'
        OTOC = np.loadtxt(path+'OTOC_Dr.dat')

        #read parameters file
        with open(path+'Parameters.dat') as f:
            params = f.read()

        #extract time step (tau) between points
        parameters = (re.findall(r'\w+', params))
        for i in range(len(parameters)):
            if parameters[i] == 'tau':
                tau = eval(parameters[i+1]+'.'+parameters[i+2])
            if parameters[i] == 'T':
                T = eval(parameters[i+1]+'.'+parameters[i+2])

        t_max,x_max = np.shape(OTOC)
        x = np.arange(0,x_max,1)
        t =np.arange(0,T+tau,tau)
        t = (len(t)/t_max)*t
        threshold = 0.1

        t_ds = np.zeros(len(x[0:x_max//2]))
        for j in x[0:x_max//2]:
            for k,time  in enumerate(t):
                if OTOC[k,j] >= threshold:
                    t_d = time
                    break
                else:
                    pass
            t_ds[j] = t_d

        ppot,pcov = opt.curve_fit(line,x[0:x_max//2],t_ds)
        sd = np.sqrt(np.diag(pcov))
        v_bs[l,m] = 1/ppot[0]

    ax.plot(taus,v_bs[:,m],label=f'$\delta j={dj}$',marker=markers[m])

fig.set_size_inches(11,6)
ax.set_xlabel('$\\tau$',size=25)
ax.set_ylabel('$v_b$',size=25)
ax.legend(fontsize=20)
ax.tick_params(labelsize=18)
plt.show()