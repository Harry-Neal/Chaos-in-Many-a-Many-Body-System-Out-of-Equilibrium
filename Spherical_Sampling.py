import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as tck

points=1000
#incorect
theta = np.pi*np.random.rand(points)
phi = 2*np.pi*np.random.rand(points)

x = np.sin(theta)*np.cos(phi)
y = np.sin(theta)*np.sin(phi)
z = np.cos(theta)

fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(2, 2, 1, projection='3d')
ax.scatter(x,y,z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('Incorect sampling')
ax.set_aspect('equal')

ax = fig.add_subplot(2, 2 , 3)
ax.scatter(phi,theta,marker='x')
ax.set_ylabel('$\\theta$')
ax.set_xlabel('$\phi$')
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))

#corect
theta = np.arccos(1-2*np.random.rand(points))
phi = 2*np.pi*np.random.rand(points)

x = np.sin(theta)*np.cos(phi)
y = np.sin(theta)*np.sin(phi)
z = np.cos(theta)

ax = fig.add_subplot(2, 2, 2, projection='3d')
ax.scatter(x,y,z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('Corect sampling')
ax.set_aspect('equal')

ax = fig.add_subplot(2, 2, 4)
ax.scatter(phi,theta,marker='x')
ax.set_ylabel('$\\theta$')
ax.set_xlabel('$\phi$')
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))

plt.show()
