
import numpy as np
import matplotlib.pyplot as plt
from data_sigmas import *

def dat (sigmas, R) :
    print('R =', R, ':')
    print('Max:', np.max(sigmas[100:]))
    print('Min:', np.min(sigmas[100:]))
    print('Avg:', np.average(sigmas[100:]))
    print('Variance:', np.var(sigmas[100:]))
    print(R, '&', np.min(sigmas[100:]), '&', np.max(sigmas[100:]), '&', np.average(sigmas[100:]), '&', np.var(sigmas[100:]))
    print()

dat(sigmas6, 6)
dat(sigmas8, 8)
dat(sigmas12, 12)
dat(sigmas20, 20)

plt.plot(range(len(sigmas20)), sigmas20, label=r'$\sigma_{R = 20}$')
plt.plot(range(len(sigmas12)), sigmas12, label=r'$\sigma_{R = 12}$')
plt.plot(range(len(sigmas8)), sigmas8, label=r'$\sigma_{R = 8}$')
plt.plot(range(len(sigmas6)), sigmas6, label=r'$\sigma_{R = 6}$')
plt.xlabel(f'$t$')
plt.ylabel(f'$\sigma$')
plt.legend()
plt.show()
