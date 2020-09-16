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
e_t = 0
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

# dt = int(input('enter dt (must be lower than dt_max)= '))
# if dt >= dt_max:
#     dt = int(input('enter dt (must be lower than dt_max)= '))
dt = 10

h = [yn_new for i in np.zeros(nxc)]
c = [sqrt(9.81 * hn) for i in np.zeros(nxc)]
q = [20 for i in np.zeros(nxc)]
b = [0.3 for i in np.zeros(nxc)]
u = [(q[i] / (b[i] * h[i])) for i, j in enumerate(np.zeros(nxc))]
r = dt / dx

time_count = int(t / dt)

cl = np.zeros([nxc, time_count])
ul = np.copy(cl)
cr = np.copy(cl)
ur = np.copy(cl)
hl = np.copy(cl)
hr = np.copy(cl)

hp = np.copy(cl)
up = np.copy(cl)

while e_t < t:
    e_t += dt
    j = int(e_t / dt)
    for i in range(1, nxc - 2):
        'left node'
        cl[i][j] = (c[i] + (r * sqrt(9.81 * h[i - 1]) * u[i] - sqrt(9.81 * h[i]) * u[i - 1])) / \
                   (1 + (r * (u[i] + sqrt(9.81 * h[i]) - (u[i - 1] + sqrt(9.81 * h[i - 1])))))
        ul[i][j] = (u[i] - r * (u[i] - u[i - 1]) * cl[i][j]) / \
                   (1 + r * (u[i] - u[i - 1]))
        hl[i][j] = h[i] - (r * (ul[i][j] + cl[i][j]) * (h[i] - h[i - 1]))

        'right node'
        cr[i][j] = (c[i] + (r * sqrt(9.81 * h[i]) * u[i + 1] - sqrt(9.81 * h[i + 1]) * u[i])) / \
                   (1 + (r * (u[i + 1] + sqrt(9.81 * h[i + 1]) - (u[i] + sqrt(9.81 * h[i])))))
        ur[i][j] = (u[i] - r * (u[i + 1] - u[i]) * cr[i][j]) / \
                   (1 + r * (u[i + 1] - u[i]))
        hr[i][j] = h[i] - (r * (ur[i][j] + cr[i][j]) * (h[i + 1] - h[i]))

    hp = [0 for i in range(1, nxc - 1)]
    up = [0 for i in range(1, nxc - 1)]

    for i in range(1, nxc - 2):
        h[0] = h[0] + (1 / 12 * 0.01)
        h[nxc - 1] = hn

        a = np.reshape(([1 / dt, cl[i][j] / (9.81 * dt), 1 / dt, cr[i][j] / (9.81 * dt)]), (2, 2))

print()
