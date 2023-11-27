import matplotlib.pyplot as plt
import numpy as np

points=1000
#incorect
theta = np.pi*np.random.rand(points)
phi = 2*np.pi*np.random.rand(points)

x = np.sin(theta)*np.cos(phi)
y = np.sin(theta)*np.sin(phi)
z = np.cos(theta)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(x,y,z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('Incorect sampling')
ax.set_aspect('equal')
plt.show()

#corect
theta = np.arccos(1-2*np.random.rand(points))
phi = 2*np.pi*np.random.rand(points)

x = np.sin(theta)*np.cos(phi)
y = np.sin(theta)*np.sin(phi)
z = np.cos(theta)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(x,y,z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('Corect sampling')
ax.set_aspect('equal')
plt.show()
