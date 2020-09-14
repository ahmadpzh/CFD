""" this is the solution for STI problem
"""

import time
import numpy as np
from numpy import sqrt
import matplotlib as plt

time1 = time.time()

## Setting
s0 = 0.001  # bed slope
b = 0.3  # flow width [m]
n = 0.03  # maning coefficient
q_l = 0  # lateral flow
q = 20 * 0.001  # water flow [m3/s]
h_g = 0.5  # head gradiant
t = 3600  # total time
L = 500  # channel length [m]
dx = 25  # space interval [m]
elasped_time = 0
epsilon = 10 ** (-5)

## Mesh Generation
nxc = int(L / dx)
nxv = nxc + 1
xv = np.zeros(nxv)
xc = np.copy(xv)
for i in range(nxc):
    xv[i + 1] = xv[i] + dx
    xc[i] = 0.5 * (xv[i + 1] + xv[i])

## Processing
yn_old = 0.1
for i in range(100):
    yn_new = ((0.03 * 0.02 * (0.3 + 2 * yn_old) ** (2 / 3)) / (0.3 ** (5 / 3) * sqrt(0.001))) ** (3 / 5)
    if yn_new - yn_old < 0.0001:
        break
    else:
        yn_old = yn_new
hn = yn_new
u_ini = q / (b * hn)
sf = ((n ** 2) * (u_ini ** 2)) / ((b + hn / (b + 2 * hn)) ** (4 / 3))
c = sqrt(9.81 * hn)
dt_max = dx / (u_ini + c)
f = c * (sf - s0)
print('dt_max is= ', dt_max)

dt = 10
h = np.zeros(nxc)
u = np.copy(h)
for i in list(range(nxc)):
    h[i] = yn_new
    c[i] = sqrt(9.81*hn)
    q[i] = 20
    b[i] = 0.3
    u[i] = q[i]/(b[i]*h[i])

print()
