import numpy as np
import re
import matplotlib.pyplot as plt
import scipy.optimize as opt

def line(x,m,c):
    return m*x + c

betas=[0.0,1.062,2.888]
djs=[0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]
v_bs=np.zeros((len(betas),len(djs)))

markers=['x','.','d']
fig = plt.figure(figsize=(11,6))
ax = plt.axes()
for p,beta in enumerate(betas):
    for m,dj in enumerate(djs):

        path = 'Figures/Fig 7/OTOC_dj_varied/'+ str(beta) +'/' + str(dj)+ '/'
        OTOC = np.loadtxt(path+'OTOC.dat')

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
        v_bs[p,m] = 1/ppot[0]
    ax.plot(djs,v_bs[p,:],label=f'$\\beta = {beta}$',linestyle='--',marker=markers[p])   

ax.set_xlabel('$\\delta J$',size=25)
ax.set_ylabel('$v_b$',size=25)
ax.legend(fontsize=20)
ax.tick_params(labelsize=17)
plt.show()