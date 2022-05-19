import numpy as np
from sympy.solvers import solve
from sympy import Symbol
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from scipy.spatial.transform import Rotation as R


def legend_without_duplicate_labels(axis):
    handles, labels = axis.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    axis.legend(*zip(*unique), loc='upper right')


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


def draw_line(v1, v2, ax, Text='', Perc=0.08, XY_plane=True, Arrow=True, Col='green', Width_Line=0.25):
    if XY_plane:
        v1 = v1.reshape(2)
        v2 = v2.reshape(2)
        if Arrow:
            v3 = np.array([0.5, 0.5])
            ax.text(*(v1 + (v2 - v1) * 0.5)-v3, Text, color=Col)
            Perc2 = Perc+0.03  # due to arrow length
            p_x1 = v1[0] + (v2[0] - v1[0]) * (0.50 - Perc)
            p_y1 = v1[1] + (v2[1] - v1[1]) * (0.50 - Perc)
            dx1  = -(v2[0] + (v1[0] - v2[0]) * (0.50 + Perc2)) + v1[0]
            dy1  = -(v2[1] + (v1[1] - v2[1]) * (0.50 + Perc2)) + v1[1]
            p_x2 = v1[0] + (v2[0] - v1[0]) * (0.50+Perc)
            p_y2 = v1[1] + (v2[1] - v1[1]) * (0.50+Perc)
            dx2  = -(v1[0] + (v2[0] - v1[0]) * (0.50+Perc2)) + v2[0]
            dy2  = -(v1[1] + (v2[1] - v1[1]) * (0.50+Perc2)) + v2[1]
            plt.arrow(p_x1, p_y1, dx1, dy1, width=Width_Line, facecolor=Col, edgecolor=Col)
            plt.arrow(p_x2, p_y2, dx2, dy2, width=Width_Line, facecolor=Col, edgecolor=Col)
        else:
            v3 = np.array([0.5, 0.5])
            ax.text(*(v1 + (v2 - v1) * 0.5) - v3, Text, color=Col)
            Perc2 = Perc  # due to arrow length
            P1 = v1 + (v2 - v1) * (0.50 - Perc)
            dP1 = -(v2 + (v1 - v2) * (0.50 + Perc2)) + v1
            P2 =  v1 + (v2 - v1) * (0.50 + Perc)
            dP2 =  -(v1 + (v2 - v1) * (0.50 + Perc2)) + v2
            ax.plot(*np.array([P1, P1+dP1]).transpose(), linewidth=10*Width_Line, color=Col)
            ax.plot(*np.array([P2, P2+dP2]).transpose(), linewidth=10*Width_Line, color=Col)
    else:
        v1 = v1.reshape(3)
        v2 = v2.reshape(3)
        v3 = np.array([0.5, 0.5, 0.5])
        ax.text(*(v1 + (v2 - v1) * 0.5) - v3, Text, color=Col)
        Perc2 = Perc  # due to arrow length
        P1 = v1 + (v2 - v1) * (0.50 - Perc)
        dP1 = -(v2 + (v1 - v2) * (0.50 + Perc2)) + v1
        P2 =  v1 + (v2 - v1) * (0.50 + Perc)
        dP2 =  -(v1 + (v2 - v1) * (0.50 + Perc2)) + v2
        ax.plot(*np.array([P1, P1+dP1]).transpose(), linewidth=10*Width_Line, color=Col)
        ax.plot(*np.array([P2, P2+dP2]).transpose(), linewidth=10*Width_Line, color=Col)


### Params
Radius         = 78.5 / 2                           # Radius of the dish
h              = Radius / 1.6                       # Height of the dish
phi            = np.radians(160)                    # Angle of the seed with respect to the x-axis
psi            = np.radians(60)                     # Angle of the pin with respect to the dish
x              = -17                                # x-coordinate of the seed
y              = -23                                # y-coordinate of the seed
margin         = 3                                  # Margin closest distance from seed to the wall
Pin_width_top  = 8                                  # Width of the pin
Pin_thick_top  = 5                                  # Thickness of the pin
Pin_width_bot  = 7                                  # Width of the pin
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
P_seed2    = P_seed - np.array([[L_lim * np.cos(phi)], [L_lim * np.sin(phi)], [0]])
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

############################## Plotting figure 1 ###########################################
minor_ticks = np.arange(-2*Radius, 2*Radius, 10)
major_ticks = np.arange(-2*Radius, 2*Radius, 5)
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
ax.scatter(*P_seed, color='k', s=50, alpha=1, label='Seed')                                # zaadje
ax.scatter(*P_1, color=col, s=scatter_size, alpha=1)                                       # pincet schaaltje onder
ax.text(P_1[0][0], P_1[1][0], P_1[2][0], '$P_1$')
ax.scatter(*P_2, color=col, s=scatter_size, alpha=1)                                       # pincet schaaltje boven
ax.scatter(*P_top, color=col, s=scatter_size, alpha=1)
ax.plot(*P_top[:, [1, 2, 0, 3, 1]], color=col, linewidth=line_width, label='Pincet frame')
ax.plot(*P_bot[:, [1, 2, 0, 3, 1]], color=col, linewidth=line_width, label='Pincet frame')
ax.text(P_2[0][0], P_2[1][0], P_2[2][0], '$P_2$')
ax.text(P_seed[0][0], P_seed[1][0]-3, P_seed[2][0], '$P_{seed} ('+str(x)+','+str(y)+')$')
for i in range(4):
    ax.text(P_top[0, i], P_top[1, i], P_top[2, i], '$P_{Top' + str(i) + '}$')
    ax.text(P_bot[0, i], P_bot[1, i], P_bot[2, i], '$P_{Bot' + str(i) + '}$')
ax.plot(*np.hstack((P_2, P_seed)), color='k', linewidth=line_width, linestyle='--')        # P2 Pseed
ax.plot(*np.hstack((P_1, P_seed)), color='k', linewidth=line_width, linestyle='--')        # P2 P1
ax.plot(*np.hstack((P_seed, P_seed2)), color='k', linewidth=3*line_width)                  # Pseed stalk
ax.plot(*np.hstack((P_seed, P_seed+np.array([[0.5*Radius], [0], [0]]))), color='k', linewidth=line_width, linestyle='--')    # Pseed line x-axis
ax.plot(8*np.cos(np.linspace(0, phi, 100)) + P_seed[0], 8*np.sin(np.linspace(0, phi, 100)) + P_seed[1], 0,
                                                                            color='k', linewidth=line_width, linestyle='--') # Pseed arc x-axis
ax.text(P_seed[0][0]+2, P_seed[1][0]+2, P_seed[2][0], '$\phi =$'+str(np.round(np.rad2deg(phi)))+'$^\circ$')                  # Pseed arc x-axis
ax.plot(*np.hstack((P_top[:, 0].reshape(3, 1), P_bot[:, 0].reshape(3, 1))), color=col, linewidth=line_width)                 # Vertical boundary
ax.plot(*np.hstack((P_top[:, 1].reshape(3, 1), P_bot[:, 1].reshape(3, 1))), color=col, linewidth=line_width)                 # Vertical boundary
ax.plot(*np.hstack((P_top[:, 2].reshape(3, 1), P_bot[:, 2].reshape(3, 1))), color=col, linewidth=line_width)                 # Vertical boundary
ax.plot(*np.hstack((P_top[:, 3].reshape(3, 1), P_bot[:, 3].reshape(3, 1))), color=col, linewidth=line_width)                 # Vertical boundary
draw_line(np.array([-Radius, 0, 0]), np.array([-Radius, 0, h]), ax, 'h', Perc=0.12,
          XY_plane=False, Arrow=False, Col='k', Width_Line=0.25)                                                         # Parameter h

if Check:
    Text = 'Seed can be grabbed :)'
else:
    Text = 'Seed cannot be grabbed :"('
ax.set_title(Text)
ax.set_xlabel('$X$-axis $[mm]$')
ax.set_ylabel('$Y$-axis $[mm]$')
ax.set_zlabel('$Z$-axis $[mm]$')
ax.set_proj_type('ortho')
legend_without_duplicate_labels(ax)
axisEqual3D(ax)

############################## Plotting figure 2 ###########################################
fig2 = plt.figure()
ax2  = plt.subplot()
ax2.quiver(np.array([0, 0]), np.array([0, 0]), np.array([1, 0]), np.array([0, 1]), scale=2.5, color='blue', label='Coordinate frame', width=0.005)
ax2.text(X[0] + U[0]*Radius, Y[0] + V[0]*Radius, '$x$', fontsize=12)
ax2.text(X[0] + U[1]*Radius, Y[0] + V[1]*Radius, '$y$', fontsize=12)                                        # Assenstelsel
ax2.plot([-Radius, Radius], [0, 0], color='k', linestyle='-', linewidth=1)                                  # Assenstelsel lines
ax2.plot([0, 0], [-Radius, Radius], color='k', linestyle='-', linewidth=1)                                  # Assenstelsel lines
ax2.plot(Radius*np.cos(np.linspace(0, 2*np.pi, 100)), Radius*np.sin(np.linspace(0, 2*np.pi, 100)),
         color='grey', linewidth=line_width, linestyle='-', label='Bakje')                                  # Bakje
ax2.plot(*np.hstack((P_seed[:2], P_seed2[:2])), color='k', linewidth=3*line_width, linestyle='-',label='Seed')  # Pseed stalk
ax2.scatter(P_seed[0], P_seed[1], color='k', s=100)                                                                 # P2 Pseed
ax2.plot(*np.hstack((P_1[:2], P_seed[:2])), color='k', linewidth=line_width, linestyle='--')                 # P2 P1
ax2.text(P_seed[0][0]-2, P_seed[1][0]-2, '$P_{seed} ('+str(x)+','+str(y)+')$')
ax2.plot([P_seed[0], P_seed[0]], [P_seed[1], 0], color='k', linestyle=':', linewidth=1)                     # P2 Pseed line x-axis
ax2.plot([P_seed[0], 0], [P_seed[1], P_seed[1]], color='k', linestyle=':', linewidth=1)                     # P2 Pseed line y-axis
ax2.plot(8*np.cos(np.linspace(0, phi, 100)) + P_seed[0], 8*np.sin(np.linspace(0, phi, 100)) + P_seed[1], 0,
        color='k', linewidth=line_width, linestyle='--')                                                    # Pseed arc x-axis
ax2.text(P_seed[0][0]+2, P_seed[1][0]+2, '$\phi =$'+str(np.round(np.rad2deg(phi)))+'$^\circ$')              # Pseed arc x-axis

draw_line(np.array([0, 0]), np.array([Radius*np.cos(np.deg2rad(135)), Radius*np.sin(np.deg2rad(135))]),
                                                        ax2, 'R', Perc=0.08, XY_plane=True, Arrow=False, Col='k', Width_Line=0.25)       # Parameter R
draw_line(np.array([0, 0]),     np.array([c, 0]),       ax2, 'c', Perc=0.08, XY_plane=True, Arrow=False, Col='green', Width_Line=0.25)   # Parameter c
draw_line(np.array([c, 0-0.5]), np.array([a+c, 0-0.5]), ax2, 'a', Perc=0.08, XY_plane=True, Arrow=False, Col='cyan', Width_Line=0.25)    # Parameter a
draw_line(np.array([c, 0]),     P_seed[:2],             ax2, 'b', Perc=0.08, XY_plane=True, Arrow=False, Col='blue', Width_Line=0.25)    # Parameter b

ax2.set_title('Sketch parameters')
ax2.set_xlabel('$X$-axis $[mm]$')
ax2.set_ylabel('$Y$-axis $[mm]$')
ax2.set_aspect('equal', 'box')
plt.minorticks_on()
plt.grid(visible=True, which='major', color='grey', linestyle='-', alpha=0.5)
plt.grid(visible=True, which='minor', color='grey', linestyle='-', alpha=0.25)

legend_without_duplicate_labels(ax2)