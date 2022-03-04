import numpy as np
from sympy.solvers import solve
from sympy import Symbol
import matplotlib.pyplot as plt


def axisEqual3D(ax):
    """
    :param ax: Generates equal axes for the 3d plot for axes ax
    :return: Plot correction
    """
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:, 1] - extents[:, 0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize / 2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)


def Draw_CoordFrame(X, Y, Z, U, V, W, ax, Text, C, L):
    ax.quiver(X, Y, Z, U, V, W, length=L, color=C)
    ax.text(X[0] + U[0] * L, Y[0] + V[0] * L, Z[0] + W[0] * L, Text + ' x', fontsize=12)
    ax.text(X[0] + U[1] * L, Y[0] + V[1] * L, Z[0] + W[1] * L, Text + ' y', fontsize=12)
    ax.text(X[0] + U[2] * L, Y[0] + V[2] * L, Z[0] + W[2] * L, Text + ' z', fontsize=12)


R     = 7.85/2
phi   = np.radians(0)
h     = R/1.6
x     = -2.1
y     = -2.3

if np.radians(-0.1) < phi < np.radians(0.1):
    Rx = (R**2-y**2)**0.5
    L = Rx-x
    print("Edge case phi = 0")
elif phi > np.radians(179.9) or phi < np.radians(-179.9):
    print("Edge case phi = 180")
    Rx = (R ** 2 - y ** 2) ** 0.5
    L  = Rx+x
else:
    a     = y/np.tan(phi)
    b     = y/np.sin(phi)
    c     = x-a
    theta = np.radians(180)-phi
    z = Symbol('x')
    z = solve(z**2 - (R**2-c**2+2*z*c*np.cos(theta)), z)
    L = max(z)-b
L_lim = h/np.tan((np.radians(60)))
print("L = ", L)
print("L>L_lim = ", L > L_lim)


############################## Plotting ###########################################
fig    = plt.figure()
ax     = fig.add_subplot(projection='3d')
theta  = np.linspace(-np.pi, 3*np.pi, 200, endpoint=True)
Height = np.array([np.zeros((len(theta))), h*np.ones((len(theta)))])
for i in range(2):                                                                  # Teken bakje
    Bakje = np.array([[R * np.cos(theta)],
                      [R * np.sin(theta)],
                      [Height[i]]])
    ax.plot(Bakje[0][0], Bakje[1][0], Bakje[2][0], color='green', linewidth=1)

U = np.array([1, 0, 0])
V = np.array([0, 1, 0])
W = np.array([0, 0, 1])
X = np.array([0, 0, 0])
Y = np.array([0, 0, 0])
Z = np.array([0, 0, 0])
Draw_CoordFrame(X, Y, Z, U, V, W, ax, '', 'blue', R)                                # Teken coordinate frame
ax.set_title('Bejo')
ax.set_xlabel('x-axis')
ax.set_ylabel('y-axis')
ax.set_zlabel('z-axis')
ax.set_proj_type('ortho')
axisEqual3D(ax)