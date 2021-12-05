
import numpy as np
import matplotlib.pyplot as plt
from data_sigmas import *

def dat (sigmas, R) :
    print('R =', R, ':')
    print('Min:', np.min(sigmas[100:])*1000)
    print('Max:', np.max(sigmas[100:])*1000)
    print('Avg:', np.average(sigmas[100:])*1000)
    print('Variance:', np.var(sigmas[100:]))
    print(R, '&', np.min(sigmas[100:]), '&', np.max(sigmas[100:]), '&', np.average(sigmas[100:]), '&', np.var(sigmas[100:]))
    print()

dat(sigmas6, 6)
dat(sigmas8, 8)
dat(sigmas12, 12)
dat(sigmas20, 20)

plt.plot(range(0, len(sigmas20), 3), sigmas20[::3], label=r'$\sigma_{R = 20}$')
plt.plot(range(0, len(sigmas12), 3), sigmas12[::3], label=r'$\sigma_{R = 12}$')
plt.plot(range(0,  len(sigmas8), 3),  sigmas8[::3], label=r'$\sigma_{R = 8}$')
plt.plot(range(0,  len(sigmas6), 3),  sigmas6[::3], label=r'$\sigma_{R = 6}$')
plt.xlabel(f'$t$')
plt.ylabel(f'$\sigma$')
plt.legend()
plt.show()
