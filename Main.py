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
    ax.quiver(X, Y, Z, U, V, W, length=L, color=C, label='Coordinate frame')
    ax.text(X[0] + U[0] * L, Y[0] + V[0] * L, Z[0] + W[0] * L, '$'+ Text + ' x$', fontsize=12)
    ax.text(X[0] + U[1] * L, Y[0] + V[1] * L, Z[0] + W[1] * L, '$'+ Text + ' y$', fontsize=12)
    ax.text(X[0] + U[2] * L, Y[0] + V[2] * L, Z[0] + W[2] * L, '$'+ Text + ' z$', fontsize=12)

def Rotz(angle):
    return R.from_euler('z', angle, degrees=False).as_matrix()



### Params
Radius         = 78.5 / 2                           # Radius of the dish
h              = Radius / 1.6                       # Height of the dish
phi            = np.radians(110)                    # Angle of the seed with respect to the x-axis
psi            = np.radians(60)                     # Angle of the pin with respect to the dish
x              = -17                                # x-coordinate of the seed
y              = -23                                # y-coordinate of the seed
margin         = 0                                  # Margin closest distance from seed to the wall
Pin_width_top  = 11                                 # Width of the pin
Pin_thick_top  = 5                                  # Thickness of the pin
Pin_width_bot  = 10                                 # Width of the pin
Pin_thick_bot  = 2                                  # Thickness of the pin
L_lim          = h/np.tan(psi)                      # Length limit to grab the seed
L_zaadje       = L_lim                              # Length seed visualisation

scatter_size   = 3
line_width     = 1

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
    # z     = (2*c*np.cos(theta)+np.sqrt(4*c**2*np.cos(theta)**2-4*(-R**2+c**2)))/2
    L     = float(max(z)-b)

### Compute box around seed
P_seed     = np.array([[x], [y], [0]])                                       # Position seed
P_1        = P_seed + np.array([[L_lim * np.cos(phi)], [L_lim * np.sin(phi)], [0]])
P_2        = P_1 + np.array([[0], [0], [h]])

### Compute box around seed top
Pin_t_proj = Pin_thick_bot / 2 * (np.cos(np.pi / 2 - psi))                               # Projected width needed for length
P_bot1     = np.array([[Pin_t_proj],  [ Pin_width_bot], [0]])
P_bot2     = np.array([[-Pin_t_proj], [-Pin_width_bot], [0]])
P_bot3     = np.array([[Pin_t_proj],  [-Pin_width_bot], [0]])
P_bot4     = np.array([[-Pin_t_proj], [ Pin_width_bot], [0]])
P_bot      = P_seed + Rotz(phi) @ np.hstack((P_bot1, P_bot2, P_bot3, P_bot4))

### Compute box around seed top
Pin_t_proj = Pin_thick_top / 2 * (np.cos(np.pi / 2 - psi))                               # Projected width needed for length
P_top1     = np.array([[Pin_t_proj],  [ Pin_width_top], [0]])
P_top2     = np.array([[-Pin_t_proj], [-Pin_width_top], [0]])
P_top3     = np.array([[Pin_t_proj],  [-Pin_width_top], [0]])
P_top4     = np.array([[-Pin_t_proj], [ Pin_width_top], [0]])
P_top      = P_2 + Rotz(phi) @ np.hstack((P_top1, P_top2, P_top3, P_top4))


P_topbot   = np.hstack((P_top, P_bot))
### Check feasibility
Ps_3_L     = np.linalg.norm(P_topbot[:2, :], axis=0)
Check      = all(Ps_3_L < (Radius-margin))
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
    ax.plot(Bakje[0][0], Bakje[1][0], Bakje[2][0], color='grey', linewidth=1, label='Bakje')

Height = np.linspace(0, h, 100)
for i in range(100):                                                                # Teken vertical bakje
    Bakje = np.array([[Radius * np.cos(theta)],
                      [Radius * np.sin(theta)],
                      [Height[i]*np.ones((len(theta)))]])
    ax.plot(Bakje[0][0], Bakje[1][0], Bakje[2][0], color='grey', linewidth=1, alpha=0.03, label='Bakje')

Rvar = np.linspace(0, Radius, 100)
for i in range(100):                                                                # Teken vertical bakje
    Bakje = np.array([[Rvar[i] * np.cos(theta)],
                      [Rvar[i] * np.sin(theta)],
                      [np.zeros((len(theta)))]])
    ax.plot(Bakje[0][0], Bakje[1][0], Bakje[2][0], color='grey', linewidth=1, alpha=0.06, label='Bakje')

U = np.array([1, 0, 0])
V = np.array([0, 1, 0])
W = np.array([0, 0, 1])
X = np.array([0, 0, 0])
Y = np.array([0, 0, 0])
Z = np.array([0, 0, 0])
Draw_CoordFrame(X, Y, Z, U, V, W, ax, '', 'blue', Radius)                                  # coordinate frame
ax.text(0, 0, 0, 'O')
ax.scatter(*P_seed, color='k', s=50, alpha=1, label='Seed')                                              # zaadje
ax.scatter(*P_1, color=col, s=scatter_size, alpha=1, label='Pincet')                                       # pincet schaaltje onder
ax.text(P_1[0][0], P_1[1][0], P_1[2][0], '$P_1$')
ax.scatter(*P_2, color=col, s=scatter_size, alpha=1, label='Pincet')                                       # pincet schaaltje boven
ax.scatter(*P_top, color=col, s=scatter_size, alpha=1, label='Pincet')
ax.plot(*P_top[:, [1, 2, 0, 3, 1]], color=col, linewidth=line_width, label='Pincet')
ax.plot(*P_bot[:, [1, 2, 0, 3, 1]], color=col, linewidth=line_width, label='Pincet')
ax.text(P_2[0][0], P_2[1][0], P_2[2][0], '$P_2$')
ax.text(P_seed[0][0], P_seed[1][0]-3, P_seed[2][0], '$P_{seed} ('+str(x)+','+str(y)+')$')
for i in range(4):
    ax.text(P_top[0, i], P_top[1, i], P_top[2, i], '$P_{Top' + str(i) + '}$')
    ax.text(P_bot[0, i], P_bot[1, i], P_bot[2, i], '$P_{Bot' + str(i) + '}$')
ax.plot(*np.hstack((P_2, P_seed)), color='k', linewidth=line_width, linestyle='--')        # P2 Pseed
ax.plot(*np.hstack((P_1, P_2)), color='k', linewidth=line_width, linestyle='--')           # P1 P2
ax.plot(*np.hstack((P_1, P_seed)), color='k', linewidth=3*line_width)                      # P1 P2
ax.plot(*np.hstack((P_seed, P_seed+np.array([[0.5*Radius], [0], [0]]))), color='k', linewidth=line_width, linestyle='--')    # Pseed line x-axis
ax.plot(8*np.cos(np.linspace(0, phi, 50)) + P_seed[0], 8*np.sin(np.linspace(0, phi, 50)) + P_seed[1], 0,
                                                                            color='k', linewidth=line_width, linestyle='--') # Pseed arc x-axis
ax.text(P_seed[0][0]+2, P_seed[1][0]+2, P_seed[2][0], '$\phi =$'+str(np.rad2deg(phi)))
ax.plot(*np.hstack((P_top[:, 0].reshape(3, 1), P_bot[:, 0].reshape(3, 1))), color=col, linewidth=line_width)                 # Vertical boundary
ax.plot(*np.hstack((P_top[:, 1].reshape(3, 1), P_bot[:, 1].reshape(3, 1))), color=col, linewidth=line_width)                 # Vertical boundary
ax.plot(*np.hstack((P_top[:, 2].reshape(3, 1), P_bot[:, 2].reshape(3, 1))), color=col, linewidth=line_width)                 # Vertical boundary
ax.plot(*np.hstack((P_top[:, 3].reshape(3, 1), P_bot[:, 3].reshape(3, 1))), color=col, linewidth=line_width)                 # Vertical boundary

if Check:
    Text = 'Seed can be grabbed :)'
else:
    Text = 'Seed cannot be grabbed :"('
ax.set_title(Text)
ax.set_xlabel('$x$-axis $[mm]$')
ax.set_ylabel('$y$-axis $[mm]$')
ax.set_zlabel('$z$-axis $[mm]$')
ax.set_proj_type('ortho')
ax.legend()
axisEqual3D(ax)

# if not (np.radians(-0.1) < phi < np.radians(0.1)) and not (phi > np.radians(179.9) or phi < np.radians(-179.9)):
#     ax.plot([0, c], [0, 0], [0, 0], color='cyan', linewidth=1)                          # Teken lengte c