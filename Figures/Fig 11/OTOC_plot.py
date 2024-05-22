import matplotlib.pyplot as plt
import numpy as np
import scipy as scpy
import scipy.optimize as opt
import re

def line(x,m,c):
    return m*x + c

path = 'OTOC_Eqlib/NNN=4/'
#read OTOC data
OTOC = np.loadtxt(path+"OTOC_Dr.dat")

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
fig1,ax1=plt.subplots(1)
fig1.set_size_inches(11,6)
print(x_max)
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
ax1.plot(x[0:x_max//2-20],tds,marker = 'x',linestyle='',label=f'$D_0={threshold}$')

ppot,pcov = opt.curve_fit(line,x[0:x_max//2-20],tds)
ax1.plot(x[0:x_max//2-20],line(x[0:x_max//2-20],*ppot),color='black')
sd = np.sqrt(np.diag(pcov))[0]
v_b = 1/ppot[0]

print(f'v_b = {v_b:.5} ± {sd:.5}')
ax1.set_xlabel('x',size=25)
ax1.set_ylabel('$t_0$',size=25)
ax1.tick_params(labelsize=18)

fig2 = plt.figure(figsize=(11,6))
ax2 = plt.axes([0.125,0.15,0.7,0.8])


#roll OTOC data along position axis so x=0 is in the middle
OTOC = np.roll(OTOC,x_max//2,axis=1)

#define new shifted x coordinates
x = np.arange(-x_max//2,x_max//2,1)

#plot shifted data
plt.pcolor(x,t,OTOC,vmin=0., vmax=1,cmap='bwr')
fig2.set_size_inches(11,6)
#plot light cone fit with 
ax2.plot(x,(1/v_b)*np.abs(x),'black',label=f"t=|j|/{v_b:.4}",linewidth=2)
ax2.set_ylim(( min(t), 200))
ax2.set_xlim((min(x),max(x)))
ax2.set_aspect('auto')
ax2.set_xlabel('j',size=30)
ax2.set_ylabel('$t$',size=30)
ax2.legend(fontsize=25,loc=3)
ax2.tick_params(labelsize=20)
#add color bar
cax = plt.axes((0.85, 0.15, 0.05, 0.8))
plt.colorbar(cax=cax)
cax.set_ylabel('D(j,t)',size=30)
cax.set_ylim(0,1)
cax.tick_params(labelsize=20)

plt.show()
