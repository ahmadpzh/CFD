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
q = 20
b = 0.3
u = [(q / (b * h[i])) for i, j in enumerate(np.zeros(nxc))]
r = dt / dx

time_count = int(t / dt)

hp = np.zeros([nxc, time_count])
up = np.zeros([nxc, time_count])

while e_t < t:
    e_t += dt
    j = int(e_t / dt)
    for i in range(1, nxc - 2):
        c_b = sqrt(9.81 * h[i])
        c_a = sqrt(9.81 * h[i - 1])
        c_c = sqrt(9.81 * h[i + 1])
        h[0] = h[0] + 0.2958

        'left node'
        cl = (c_b + (r * (c_a * u[i] - c_b * u[i - 1]))) / (1 + (r * (u[i] + c_b - (u[i - 1] + c_a))))
        ul = (u[i] - r * (u[i] - u[i - 1]) * cl) / (1 + r * (u[i] - u[i - 1]))
        hl = h[i] - (r * (ul + cl) * (h[i] - h[i - 1]))

        'right node'
        cr = (c_b + (r * (c_b * u[i + 1] - c_c * u[i]))) / (1 + (r * (u[i + 1] - c_c - (u[i] - c_b))))
        ur = (u[i] + r * (u[i + 1] - u[i]) * cr) / (1 + r * (u[i + 1] - u[i]))
        hr = h[i] - (r * (ur - cr) * (h[i + 1] - h[i]))

    hp = [0 for i in range(1, nxc - 1)]
    up = [0 for i in range(1, nxc - 1)]

    for i in range(1, nxc - 2):
        h[0] = h[0] + (1 / 12 * 0.01)
        h[nxc - 1] = hn

        a = np.reshape((1 / dt, cl / (9.81 * dt), 1 / dt, cr / (9.81 * dt)), (2, 2))
        b = np.reshape((hl / dt + (cl * ul / (9.81 * dt)) + (
                    c_b * (((n ** 2 * ul ** 2) / (((b * h[i]) / (b + 2 * h[i])) ** (4 / 3))) - s0)),
                       hr / dt + (cr * ur / (9.81 * dt)) + (
                                   c_b * (s0 - ((n ** 2 * ur ** 2) / (((b * h[i]) / (b + 2 * h[i])) ** (4 / 3)))))),
                       (2, 1))

        c = np.linalg.solve(a, b)

        hp[i] = c[0]
        up[i] = c[1]
print()
print()
