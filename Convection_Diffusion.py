"""this program solves 2D Convection-Diffusion equation"""

import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

time1 = time.time()

# # Setting
l = 10  # length
h = l / 5  # height
n_x = 10  # grid in x-axis
n_y = 10  # grid in y-axis
n_x_max = n_x + 1
n_y_max = n_y + 1
x0 = 0
y0 = 0
dx = (l / n_x)
dy = (h / n_y)
ro = 1
u_mean = 1
u_max = 1.5 * u_mean
reynolds = 10
pr = 1
mu = ro * u_mean * h / reynolds
k_cp = mu / pr
# variables (Acceleration Part)
x = np.reshape(np.zeros(n_y_max * n_x_max), [n_y_max, n_x_max])
y = np.copy(x)
gamma = np.copy(x)
gamma_e = np.copy(x)
gamma_w = np.copy(x)
gamma_n = np.copy(x)
gamma_s = np.copy(x)
u = np.copy()
v = np.copy()

for i in range(n_y_max):
    for j in range(n_x_max):
        x[i][j] = float(j) * dx
        y[i][j] = float(i) * dy
        gamma[i][j] = k_cp
        u[i][j] = u_max*(1 - (2*y/))
        v[i][j] = 0

d_a = 100
omega = 0.7
error = 1000
std_error = 0.001
iteration = 0

# # Initial Condition
# 50[m] for every cells except top boundary
T_old = np.reshape([float(d_a) for i in np.zeros(n_y_max * n_x_max)], [n_y_max, n_x_max])
# creating zero matrix for process acceleration
T_new = np.reshape([float(d_a) for i in np.zeros(n_y_max * n_x_max)], [n_y_max, n_x_max])

# # Boundary Condition
# top boundary
T_old[0] = 100  # Dirichlet boundary
T_new[0] = 100  # Dirichlet boundary

# left boundary
q = 0.1  # newmann boundary [m/day]

# # Processing
while error > std_error:
    iteration += 1

    for i in range(1, n_y_max - 1):
        for j in range(1, n_x_max - 1):
            gamma_e = (gamma[i][j + 1] + gamma[i][j]) / 2
            gamma_w = (gamma[i][j - 1] + gamma[i][j]) / 2
            gamma_n = (gamma[i - 1][j] + gamma[i][j]) / 2
            gamma_s = (gamma[i + 1][j] + gamma[i][j]) / 2
            dx_e = x[i][j + 1] - x[i][j]
            dx_w = x[i][j] - x[i][j - 1]
            dy_n = y[i][j] - y[i - 1][j]
            dy_s = y[i + 1][j] - y[i][j]

            area_e = y[i + 1][j + 1] - y[i][j + 1]
            area_w = y[i + 1][j] - y[i][j]
            area_n = x[i + 1][j + 1] - x[i + 1][j]
            area_s = x[i][j + 1] - x[i][j]

            a_known_e = gamma_e * area_e / dx_e
            a_known_w = gamma_w * area_w / dx_w
            a_known_n = gamma_n * area_n / dy_n
            a_known_s = gamma_s * area_s / dy_s
            s_p = 0.
            s_u1, s_u2, s_u3, s_u4 = 0., 0., 0., 0.

            # top boundary
            if i == 1:
                a_known_n = 2 * a_known_n
                # s_u2 = a_known_n * T_old[i + 1][j]
                # s_p = -a_known_n
                s_u2 = 0
                s_p = 0
                a_known_n = 0

            # left_boundary (constant_flow)
            if j == 1:
                a_known_w = 0
                s_u3 = q * d_a * dy

            # right boundary (no flow)
            if j == n_x - 2:
                a_known_e = 0
                s_u4 = 0

            # bottom boundary (no flow)
            if i == n_y_max - 2:
                a_known_s = 0
                s_u1 = 0

            source = s_u1 + s_u2 + s_u3 + s_u4 - s_p
            a_p = a_known_e + a_known_w + a_known_n + a_known_s - s_p

            T_new[i][j] = omega * (
                    a_known_e * T_old[i][j + 1] + a_known_w * T_old[i][j - 1] + a_known_n * T_old[i + 1][j] +
                    a_known_s * T_old[i - 1][j] + source) / a_p + (1 - omega) * T_old[i][j]

    for j in range(n_x_max):
        T_new[-1][j] = T_new[-2][j]

    for i in range(n_y_max):
        T_new[i][-1] = T_new[i][-2]
        T_new[i][0] = T_new[i][1] + q

    error = np.linalg.norm(T_new - T_old, 2)
    print('\nL2Norm = %0.10f' % error)

    print('h center= %0.10f' % (T_new[5][5] - T_old[5][5]))

    print('iteration = ', iteration)

    # plt.contourf(T_new)
    # plt.gca().invert_yaxis()
    # plt.axis('off')
    # plt.grid()
    # plt.colorbar().ax.set_ylabel('[m]')
    # plt.pause(0.0001)
    # plt.show(block=False)
    # plt.clf()

    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # surf = ax.plot_surface(x, y, T_new, cmap=cm.viridis)
    # fig.colorbar(surf)
    # plt.show(block=False)
    # plt.pause(1)
    # plt.close()
    for i in range(n_y_max):
        for j in range(n_x_max):
            T_old[i][j] = T_new[i][j]

print('L2Norm = ', error)
print('iteration = ', iteration)

plt.contourf(T_new)
plt.gca().invert_yaxis()
plt.colorbar().ax.set_ylabel('[m]', rotation=270)
fig = plt.figure()
ax = fig.gca(projection='3d')
fig.gca().invert_yaxis()
surf = ax.plot_surface(x, y, T_new, cmap=cm.viridis)
fig.colorbar(surf)
# plt.savefig('Final_Result.png')
plt.show(block=False)

time2 = time.time()

print('\nTotal time = ', (time2 - time1) / 60)

print()
