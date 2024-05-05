import numpy as np
import matplotlib.pyplot as plt
import re
import scipy.optimize as opt

def line(x,m,c):
    return m*x + c


path = 'OTOC_disorder/'

djs= [0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]
vbs = np.zeros(len(djs))

for m,dj in enumerate(djs):

    path = 'OTOC_disorder/' + str(dj) +'/'
    OTOC = np.loadtxt(path+'OTOC.dat')

    #read parameters file
    with open(path+'Parameters.dat') as f:
        params = f.read()

    #extract time step (tau) between points
    parameters = (re.findall(r'\w+', params))
    for i in range(len(parameters)):
        if parameters[i] == 'dt':
            dt = eval(parameters[i+1]+'.'+parameters[i+2])
        if parameters[i] == 'T':
            T = eval(parameters[i+1]+'.'+parameters[i+2])

    t_max,x_max = np.shape(OTOC)

    x = np.arange(0,x_max,1)
    t = np.arange(0,T+dt,dt)

    N=100
    thresholds = np.linspace(0.1,0.5,N)
    v_b = np.zeros(len(thresholds))

    for i,threshold in enumerate(thresholds):
        t_ds = np.zeros(len(x[0:x_max//2]))
        for j in x[0:x_max//2]:
            for k,time in enumerate(t):
                if OTOC[k,j] >= threshold:
                    t_d = time
                    break
                else:
                    pass
            t_ds[j] = t_d

        ppot,pcov = opt.curve_fit(line,x[0:x_max//2],t_ds)
        sd = np.sqrt(np.diag(pcov))
        v_b[i] = 1/ppot[0]
    vbs[m] = np.average(v_b)

fig,ax = plt.subplots(1)
ax.plot(djs,vbs,marker='x')
ax.set_xlabel('$\delta j$',size = 30)
ax.set_ylabel('$v_b$', size =30)
ax.tick_params(labelsize=20)

plt.show()