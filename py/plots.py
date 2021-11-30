
import numpy as np
import matplotlib.pyplot as plt
import eigen_v


#
# Utilities to plot results from simulations
#



# Increase font size on figures
smallFig = True
plt.rcParams.update({'font.size': 12 if smallFig else 14})


# Import simulation results from data.py
from data import y, u_x, t, omega


# Compute analytical results in discrete steps
y[0] += 0.5
y[len(y)-1] -= 0.5
pulsatile = omega > -1
a = y[len(y) - 1]

y_lin = np.linspace(-a, a, len(y) * 10)

# Constant-force analytical solutions
u_x_a = eigen_v.u_x(y_lin, t, a=a, nu=0.5)
u_x_inf = 2e-5 - 5e-8 * y_lin ** 2

# Pulsatile analytical solution
u_tilde = eigen_v.u_x_tilde(y_lin, t, a=a, omega=omega, nu=0.5)

# Plot data
plt.figure('t = ' + str(t), figsize=(5.75, 6) if smallFig else (10, 8))
if pulsatile:
    plt.plot(u_tilde, y_lin, '-', label=r'$\tilde{u}$')
else:
    plt.plot(u_x_a / 2e-5, y_lin, '-', label=r'$u_x$')
    plt.plot(u_x_inf / 2e-5, y_lin, 'k--', label=r'$u_x^{\infty}$', alpha=0.5)
plt.plot(u_x / (1 if pulsatile else 2e-5), y, 'r+', label=r'$u_x^{LBM}$')
plt.legend()
plt.ylabel(r'$y$ (lattice constants)')
if pulsatile:
    plt.xlabel(r'velocity')
else:
    plt.xlabel(r'$\frac{velocity}{2 \times 10^{-5}}$')


# Auto-generated max velocities at various time-steps
if False:
    velocities = [0, 0.0000024463826215863682, 0.000004835599960464384, 0.0000069738458324549756, 0.000008834684428318178, 0.00001043674051100514, 0.00001181159719092598, 0.000012991627949987935, 0.000014003945489519874, 0.000014871920512076276, 0.000015616296381910727, 0.00001625476306419426, 0.00001680230560784447, 0.00001727185879853915, 0.000017674560424575315, 0.000018019926543047403, 0.00001831611149538636, 0.00001857012114531806, 0.000018787963355521025, 0.000018974786704345973, 0.000019135007531733296, 0.000019272414322325785, 0.00001939025565494477, 0.000019491317323770257, 0.000019577988628323475, 0.00001965231868219845, 0.00001971606476303125, 0.000019770733928661652, 0.000019817618662074153, 0.000019857827403553654, 0.000019892310755961723, 0.000019921883966456553, 0.0000199472461994108, 0.000019968997062738663, 0.000019987650784877427, 0.000020003648372995665, 0.000020017368038383, 0.000020029134138201497, 0.00002003922484380626, 0.00002004787871727025, 0.00002005530035157387]
    labels = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000]
    t_vals = np.linspace(0, labels[len(labels)-1], 500)
    u_x_vals = eigen_v.u_x(0, t_vals)
    plt.figure('max_vel_by_timestep', figsize=(10, 8))
    plt.plot(t_vals, u_x_vals, '-', label=r'$u_x$')
    plt.plot([0, labels[len(labels) - 1]], [2e-5, 2e-5], 'k--', label=r'$u_x^{\infty}(0)$', alpha=0.5)
    plt.plot(labels, velocities, 'r+', label=r'$u_x^{LBM}(0, t)$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'Velocity at $y=0$')
    plt.legend()


plt.show()
