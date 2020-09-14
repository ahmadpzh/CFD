""" this is the solution for STI problem
"""

import time
import numpy as np
import matplotlib as plt

time1 = time.time()

## Setting
s = 0.001               # bed slope
b = 0.3                 # flow width [m]
n = 0.03                # maning coefficient
q_l = 0                 # lateral flow
q = 20 * 0.001          # water flow [m3/s]
h_g = 0.5               # head gradiant
t = 3600                # total time
l = 500                 # channel length [m]
dx = 25                 # space interval [m]
elasped_time = 0
epsilon = 10 ** (-5)

## Mesh Generation
nxc = l/dx
nxv = nxc + 1
xv = np.zeros(nxv)
for i in range(nxc):
    xv[i+1] = xv[i] + dx