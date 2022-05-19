import numpy as np

### Params
Radius         = 85 / 2                           # Radius of the dish
h              = 40                       # Height of the dish
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
    # z     = [np.sqrt(Radius**2-c**2*np.sin(theta)**2)+c*np.cos(theta),
    #          -np.sqrt(Radius**2-c**2*np.sin(theta)**2)+c*np.cos(theta)]
    # L     = float(max(z)-b)

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
Rotz       = np.array([[np.cos(phi), -np.sin(phi), 0],
                       [np.sin(phi), np.cos(phi), 0],
                       [0, 0, 1]])
P_bot      = P_seed + Rotz @ np.hstack((P_bot1, P_bot2, P_bot3, P_bot4))

### Compute box around seed top
Pin_t_proj = Pin_thick_top / 2 * (np.cos(np.pi / 2 - psi))                               # Projected width needed for length
P_top1     = np.array([[Pin_t_proj],  [ Pin_width_top], [0]])
P_top2     = np.array([[-Pin_t_proj], [-Pin_width_top], [0]])
P_top3     = np.array([[Pin_t_proj],  [-Pin_width_top], [0]])
P_top4     = np.array([[-Pin_t_proj], [ Pin_width_top], [0]])
P_top      = P_2 + Rotz @ np.hstack((P_top1, P_top2, P_top3, P_top4))


P_topbot   = np.hstack((P_top, P_bot))
### Check feasibility
Ps_3_L     = np.linalg.norm(P_topbot[:2, :], axis=0)
Check      = all(Ps_3_L < (Radius-margin))
if Check:
    print('Seed can be grabbed :)')
else:
    print('Seed cannot be grabbed :"(')