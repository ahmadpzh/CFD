"""this program solves 2D equation"""

import time
import numpy as np
import matplotlib.pyplot as plt

# time1 = time.time()

# # Setting
l_x = 600  # length
l_y = 400  # width
std_error = 0.001
n_x = 30  # grid in x-axis
n_y = 20  # grid in y-axis
n_x_max = n_x+1
n_y_max = n_y+1
x0 = 0
y0 = 0
dx = (l_x / n_x_max)
dy = (l_y / n_y_max)
x = np.reshape(np.zeros(n_y_max * n_x_max), [n_y_max, n_x_max])
y = np.copy(x)
h = np.copy(x)

for i in range(n_y_max):
    for j in range(n_x_max):
        x[i][j] = float(j) * dx
        y[i][j] = float(i) * dy
depth_average = 50
h[0] = 100









print()
