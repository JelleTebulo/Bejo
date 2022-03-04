import numpy as np
from sympy.solvers import solve
from sympy import Symbol
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R


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

def Rotz(angle):
    return R.from_euler('z', angle, degrees=False).as_matrix()



### Params
Radius     = 7.85 / 2                           # Radius of the dish
h          = Radius / 1.6                       # Height of the dish
phi        = np.radians(120)                    # Angle of the seed with respect to the x-axis
psi        = np.radians(60)                     # Angle of the pin with respect to the dish
x          = -1.7                               # x-coordinate of the seed
y          = -2.3                               # y-coordinate of the seed
Pin_width  = 1                                  # Width of the pin
Pin_thick  = 0.3                                # Thickness of the pin
L_lim      = h/np.tan((psi))                    # Length limit to grab the seed
L_zaadje   = L_lim                              # Length seed visualisation

scatter_size = 3
line_width   = 1

### Compute L
if np.radians(-0.1) < phi < np.radians(0.1):
    Rx = (Radius ** 2 - y ** 2) ** 0.5
    L = Rx-x
    print("Edge case phi = 0")
elif phi > np.radians(179.9) or phi < np.radians(-179.9):
    print("Edge case phi = 180")
    Rx = (Radius ** 2 - y ** 2) ** 0.5
    L  = Rx+x
else:
    a     = y/np.tan(phi)
    b     = y/np.sin(phi)
    c     = x-a
    theta = np.radians(180)-phi
    z     = Symbol('x')
    z     = solve(z ** 2 - (Radius ** 2 - c ** 2 + 2 * z * c * np.cos(theta)), z)
    L     = float(max(z)-b)

### Compute box around seed
Pin_t_proj = Pin_thick/2*(np.cos(np.pi/2-psi))                               # Projected width needed for length
P_seed     = np.array([[x], [y], [0]])                                       # Position seed
P_1        = P_seed + np.array([[L_lim * np.cos(phi)], [L_lim * np.sin(phi)], [0]])
P_2        = P_1 + np.array([[0], [0], [h]])
P_31       = np.array([[Pin_t_proj], [Pin_width], [0]])
P_32       = np.array([[-Pin_t_proj], [-Pin_width], [0]])
P_33       = np.array([[Pin_t_proj],  [-Pin_width], [0]])
P_34       = np.array([[-Pin_t_proj], [Pin_width],  [0]])
Ps_3       = P_2 + Rotz(phi) @ np.hstack((P_31, P_32, P_33, P_34))
Ps_4       = P_seed + Rotz(phi) @ np.hstack((P_31, P_32, P_33, P_34))
Ps         = np.hstack((Ps_3, Ps_4))
### Check feasibility
Ps_3_L     = np.linalg.norm(Ps[:2, :], axis=0)
Check      = all(Ps_3_L < Radius)
if Check:
    col = 'green'
else:
    col = 'red'

############################## Plotting ###########################################
fig    = plt.figure()
ax     = fig.add_subplot(projection='3d')
theta  = np.linspace(-np.pi, 3*np.pi, 200, endpoint=True)
Height = np.array([np.zeros((len(theta))), h*np.ones((len(theta)))])
for i in range(2):                                                                  # Teken bakje 2 lijnen
    Bakje = np.array([[Radius * np.cos(theta)],
                      [Radius * np.sin(theta)],
                      [Height[i]]])
    ax.plot(Bakje[0][0], Bakje[1][0], Bakje[2][0], color='grey', linewidth=1)

Height = np.linspace(0, h, 100)
for i in range(100):                                                                # Teken vertical bakje
    Bakje = np.array([[Radius * np.cos(theta)],
                      [Radius * np.sin(theta)],
                      [Height[i]*np.ones((len(theta)))]])
    ax.plot(Bakje[0][0], Bakje[1][0], Bakje[2][0], color='grey', linewidth=1, alpha=0.02)

Rvar = np.linspace(0, Radius, 100)
for i in range(100):                                                                # Teken vertical bakje
    Bakje = np.array([[Rvar[i] * np.cos(theta)],
                      [Rvar[i] * np.sin(theta)],
                      [np.zeros((len(theta)))]])
    ax.plot(Bakje[0][0], Bakje[1][0], Bakje[2][0], color='grey', linewidth=1, alpha=0.04)

U = np.array([1, 0, 0])
V = np.array([0, 1, 0])
W = np.array([0, 0, 1])
X = np.array([0, 0, 0])
Y = np.array([0, 0, 0])
Z = np.array([0, 0, 0])
Draw_CoordFrame(X, Y, Z, U, V, W, ax, '', 'blue', Radius)            # Teken coordinate frame
ax.text(0, 0, 0, 'O')
ax.scatter(*P_seed, color='cyan', s=30, alpha=1)                     # Teken zaadje
ax.scatter(*P_1, color=col, s=scatter_size, alpha=1)                 # Teken pincet schaaltje onder
ax.text(P_1[0][0], P_1[1][0], P_1[2][0], '$P_1$')
ax.scatter(*P_2, color=col, s=scatter_size, alpha=1)                 # Teken pincet schaaltje boven
ax.scatter(*Ps_3, color=col, s=scatter_size, alpha=1)
ax.plot(*Ps_3[:, [1, 2, 0, 3, 1]], color=col, linewidth=line_width)
ax.plot(*Ps_4[:, [1, 2, 0, 3, 1]], color=col, linewidth=line_width)
ax.text(P_2[0][0], P_2[1][0], P_2[2][0], '$P_2$')
ax.text(P_seed[0][0], P_seed[1][0], P_seed[2][0], '$(P_{seed}'+str(x)+','+str(y)+')$')
for i in range(4):
    ax.text(Ps_3[0, i], Ps_3[1, i], Ps_3[2, i], '$P_{3'+str(i)+'}$')
    ax.text(Ps_4[0, i], Ps_4[1, i], Ps_4[2, i], '$P_{4'+str(i)+'}$')
ax.plot(*np.hstack((P_2, P_seed)), color='k', linewidth=line_width)
ax.plot(*np.hstack((P_1, P_2)), color='k', linewidth=line_width)
ax.plot(*np.hstack((P_1, P_seed)), color='k', linewidth=line_width)
ax.plot(*np.hstack((Ps_3[:, 0].reshape(3, 1), Ps_4[:, 0].reshape(3, 1))), color=col, linewidth=line_width)
ax.plot(*np.hstack((Ps_3[:, 1].reshape(3, 1), Ps_4[:, 1].reshape(3, 1))), color=col, linewidth=line_width)
ax.plot(*np.hstack((Ps_3[:, 2].reshape(3, 1), Ps_4[:, 2].reshape(3, 1))), color=col, linewidth=line_width)
ax.plot(*np.hstack((Ps_3[:, 3].reshape(3, 1), Ps_4[:, 3].reshape(3, 1))), color=col, linewidth=line_width)

if Check:
    Text = 'Seed can be grabbed :)'
else:
    Text = 'Seed cannot be grabbed :('
ax.set_title(Text)
ax.set_xlabel('$x$-axis')
ax.set_ylabel('$y$-axis')
ax.set_zlabel('$z$-axis')
ax.set_proj_type('ortho')
axisEqual3D(ax)

# if not (np.radians(-0.1) < phi < np.radians(0.1)) and not (phi > np.radians(179.9) or phi < np.radians(-179.9)):
#     ax.plot([0, c], [0, 0], [0, 0], color='cyan', linewidth=1)                          # Teken lengte c