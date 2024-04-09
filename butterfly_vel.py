import numpy as np
import re
import matplotlib.pyplot as plt
import scipy.optimize as opt

def line(x,m,c):
    return m*x + c
taus = np.arange(0.5,5.5,0.5)
for tau_cur in taus:
    path = 'OTOCs/' + str(tau_cur) + '/'
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

    t_max,x_max = np.shape(OTOC)
    x = np.arange(0,x_max,1)
    t = np.arange(0,T,tau)
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

    print(f'tau = {tau_cur}: v_b = {np.average(v_b):.5} ± {np.std(v_b):.5}')