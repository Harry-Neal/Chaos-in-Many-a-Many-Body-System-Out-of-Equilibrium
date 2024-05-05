import numpy as np
import re
import matplotlib.pyplot as plt
import scipy.optimize as opt

def line(x,m,c):
    return m*x + c

taus = [0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0]
djs=[0.05,0.1,0.25,0.35,0.5]
v_bs=np.zeros((len(taus),len(djs)))
fig = plt.figure(figsize=(10,8))
ax = plt.axes()

for m,dj in enumerate(djs):
    for l,tau_cur in enumerate(taus):
        path = 'OTOCsdj'+ str(dj) +'/' + str(tau_cur) + '/'
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
        t=(len(t)/t_max)*t
        N=100
        thresholds = np.linspace(0.1,0.5,N)
        v_b = np.zeros(len(thresholds))

        for i,threshold in enumerate(thresholds):
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
            v_b[i] = 1/ppot[0]
        v_bs[l,m] = np.average(v_b)
    
    ax.plot(taus,v_bs[:,m],label=f'dj={dj}',marker='x')

ax.set_xlabel('$\\tau$',size=30)
ax.set_ylabel('$v_b$',size=30)
ax.legend(fontsize=20)
ax.tick_params(labelsize=20)
plt.show()