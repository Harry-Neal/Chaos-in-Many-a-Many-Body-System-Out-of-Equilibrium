import numpy as np
data = np.loadtxt("OTOC.dat")
t = data[:,0]
x = data[:,1]
OTOC = data[:,2]
runs = 100
ssize = int(x.max() - x.min() + 1)
t_range = (t.max() - t.min() + 1)
dt = t[ssize] - t[0]
len_ts = int(t_range/dt)

runs_avg = np.zeros(ssize*len_ts)
for i in range(runs):
    runs_avg = runs_avg + OTOC[i*ssize*len_ts : ssize*len_ts + i*ssize*len_ts]

runs_avg = runs_avg/runs

OTOC_avg = np.zeros((len_ts,ssize))
for i in range(len_ts):
    for j in range(ssize):
        OTOC_avg[i,j] = runs_avg[j + i*ssize]

np.savetxt('OTOC_avg.dat',OTOC_avg,fmt='%.4f')
params = np.array([runs,dt])
np.savetxt('plotting_params_OTOC.dat',params,fmt='%i')