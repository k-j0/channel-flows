
import numpy as np
import matplotlib.pyplot as plt


#
# Utilities to compute analytical results for v and u_x in channel, to be compared with numerical results
#


def v(y: float, t: float, f_x: float = 5e-8, a: float = 20, nu: float = 0.5, iterations: int = 1000):
    """
        Computes v(y, t) using eigenfunction expansion
        Parameters:
            y (float): Y-coordinate within the channel, in -a..a
            t (float): Time at which to evaluate v
            f_x (float): Constant force, applied in the direction of the flow (towards x)
            a (float): Half-width of the channel
            nu (float): Kinematic viscosity
            iterations (int): Where to stop the infinite eigenfunction expansion series
        Returns:
            (float): Approximation to v(y, t) = u_x(y, t) - u_x^\infty(y, t)
    """

    sum = 0
    for j in range(1, iterations+1):
        if j%2 == 0:
            continue
        xi_j = np.sin(j * np.pi * (y + a)/(2.0 * a)) / np.sqrt(a)
        lambda_j = -nu * (j * np.pi / (2.0 * a)) ** 2
        t_term = np.exp(lambda_j * t)
        v_j = -f_x * ((8*a**3)/(j**3*np.pi**3) - (4*a**3 * (np.sin(j * np.pi) + (2 * (-1)**j)/(j * np.pi)))/(j**2*np.pi**2)) / (nu * np.sqrt(a))
        sum += v_j * t_term * xi_j
    return sum

def u_x(y: float, t: float, f_x: float = 5e-8, a: float = 20, nu: float = 0.5, iterations: int = 1000):
    """
        Computes u_x(y, t) using eigenfunction expansion of v = u_x - u_x^inf
        Parameters:
            y (float): Y-coordinate within the channel, in -a..a
            t (float): Time at which to evaluate v
            f_x (float): Constant force, applied in the direction of the flow (towards x)
            a (float): Half-width of the channel
            nu (float): Kinematic viscosity
            iterations (int): Where to stop the infinite eigenfunction expansion series
        Returns:
            (float): Approximation to u_x as u_x^\infty(y, t) + v(y, t)
    """

    u_x_inf = f_x * (a ** 2 - y ** 2) / (2.0 * nu)
    return u_x_inf + v(y, t, f_x, a, nu, iterations)


def main():
    t_vals = np.linspace(0, 3000, 500)
    u_x_vals = u_x(0, t_vals)
    
    plt.rcParams.update({'font.size': 14})

    # plt.figure()
    # plt.plot(t_vals, u_x_vals)
    # plt.ylabel('u_x(0, t)')
    # plt.xlabel('t')

    plt.figure('v_0_0-func-of-j_max', figsize=(10, 8))
    iteration_vals = range(5, 99, 2)
    v_vals = [0] * len(iteration_vals)
    for i, iter in enumerate(iteration_vals):
        v_vals[i] = v(0, 0, iterations=iter)
        if i == len(iteration_vals) - 1:
            print('v_(', iter, ') + uxinf =', v_vals[i] + 2e-5)
    plt.plot(iteration_vals, v_vals, '.-')
    plt.ylabel(r'$v(0, 0)$')
    plt.xlabel(r'$j_{max}$')

    plt.show()

if __name__ == '__main__':
    main()
