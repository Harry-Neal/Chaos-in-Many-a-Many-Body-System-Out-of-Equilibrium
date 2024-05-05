import matplotlib.pyplot as plt
import numpy as np
import scipy as scpy
import re

path = 'Next Nearest Neighbour'
#read OTOC data
OTOC = np.loadtxt(path+"/OTOC_Dr.dat")

#read parameters file
with open(path+'/Parameters.dat') as f:
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

#declare range for OTOC data
t_max,x_max = np.shape(OTOC)
x = np.arange(0,x_max,1)
t = np.arange(0,T+tau,tau)

#=============FIRST FIGURE===============
#=========CAN IGNORE FOR NOW============

# fig = plt.figure()
# ax2 = plt.axes([0.1,0.1,0.8,0.8])
# ax1=plt.axes([0.15,0.15,0.4,0.4])
# #Plot of the propoagtion of the chaotic wavefront from x=0

# #index at selected times
# times = np.array([40,50,60,70,80,90,100])

# #index OTOC at each time step and plot as function of position
# for time in times:
#    t_sample = int(np.argwhere(t == time)[0])
#    ax1.plot(x[0:200],OTOC[t_sample,0:200],label=f't={t_sample}')

# ax1.set_xlabel('x')
# ax1.set_ylabel('D(x,t)')

# #rescale data to universal curve 
# y=np.log(OTOC[:,1:200]/(epsilon**2))
# #index recaled OTOC data at each time step and plot as function of position
# for time in times:
#    t_sample = int(np.argwhere(t == time)) 
#    ax2.plot(x[1:200]/time , y[t_sample,:]/(2*time),label=f't={time}')

# #plot rescaled data
# ax2.set_xlabel('x/t')
# ax2.set_ylabel('$ln(D(x,t)/ \epsilon^2 )/2t$')
# ax2.legend()

#============SECOND FIGURE===============
fig = plt.figure(figsize=(10,8))
ax3 = plt.axes([0.1,0.1,0.7,0.8])


#roll OTOC data along position axis so x=0 is in the middle
OTOC = np.roll(OTOC,x_max//2,axis=1)

#define new shifted x coordinates
x = np.arange(-x_max//2,x_max//2,1)

#plot shifted data
plt.pcolor(x,t,OTOC,vmin=0., vmax=1,cmap='bwr')
#set aspect and labels


#plot light cone fit with 
v_b = 3.7715
ax3.plot(x,(1/v_b)*np.abs(x) + 5,'black',label=f"t=|x|/{v_b}",linewidth=2)
ax3.set_ylim(( min(t), max(t)))
ax3.set_xlim(-250,250)
ax3.set_aspect('auto')
ax3.set_xlabel('x',size=30)
ax3.set_ylabel('$t$',size=30)
ax3.legend(fontsize=20)
ax3.tick_params(labelsize=20)
#add color bar
cax = plt.axes((0.85, 0.1, 0.05, 0.8))
plt.colorbar(cax=cax)
cax.set_ylabel('D(x,t)',size=30)
cax.set_ylim(0,1)
cax.tick_params(labelsize=20)

plt.show()
