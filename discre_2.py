"""this program solves 2D equation"""

import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

time1 = time.time()

# # Setting
l_x = 600  # length
l_y = 400  # height
n_x = 10  # grid in x-axis
n_y = 10  # grid in y-axis
n_x_max = n_x + 1
n_y_max = n_y + 1
x0 = 0
y0 = 0
dx = (l_x / n_x)
dy = (l_y / n_y)
# variables (Acceleration Part)
x = np.reshape(np.zeros(n_y_max * n_x_max), [n_y_max, n_x_max])
y = np.copy(x)
gamma = np.copy(x)
gamma_e = np.copy(x)
gamma_w = np.copy(x)
gamma_n = np.copy(x)
gamma_s = np.copy(x)

for i in range(n_y_max):
    for j in range(n_x_max):
        x[i][j] = float(j) * dx
        y[i][j] = float(i) * dy
        gamma[i][j] = 1000

d_a = 50
omega = 1
error = 1000
std_error = 0.001
iteration = 0

# # Initial Condition
# 50[m] for every cells except top boundary
h_old = np.reshape([float(d_a) for i in np.zeros(n_y_max * n_x_max)], [n_y_max, n_x_max])
# creating zero matrix for process acceleration
h_new = np.reshape([float(d_a) for i in np.zeros(n_y_max * n_x_max)], [n_y_max, n_x_max])

# # Boundary Condition
# top boundary
h_old[0] = 100  # Dirichlet boundary
h_new[0] = 100  # Dirichlet boundary

# left boundary
q = 0.1  # newmann boundary [m/day]

# # Governing Equation (Water Head Distribution)
# def whd(a_p=a_known_e + a_known_w + a_known_n + a_known_s, side_coefficient_w=a_known_w, side_coefficient_e=a_known_e,
#         side_coefficient_n=a_known_n, side_coefficient_s=a_known_s, source=0.):
#     h_p = ((1 / a_p) * (side_coefficient_w * h_old[i][j - 1] + side_coefficient_e * h_old[i][j + 1] +
#                         side_coefficient_n * h_old[i - 1][j] + side_coefficient_s * h_old[i + 1][j] + source))
#     return h_p


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
            # area_e = dy_n * d_a
            # area_w = dy_n * d_a
            # area_n = dx_e * d_a
            # area_s = dx_w * d_a
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
            if i == n_y_max - 2:
                a_known_n = 2 * a_known_n
                s_u1 = a_known_n * h_old[i + 1][j]
                s_p = -a_known_n
                a_known_n = 0

            # left_boundary (constant_flow)
            if j == 1:
                a_known_w = 0
                s_u2 = q * d_a * dy

            # right boundary (no flow)
            if j == n_x - 2:
                a_known_e = 0
                s_u3 = 0

            # bottom boundary (no flow)
            if i == 1:
                a_known_s = 0
                s_u4 = 0

            source = s_u1 + s_u2 + s_u3 + s_u4 - s_p
            a_p = a_known_e + a_known_w + a_known_n + a_known_s - s_p

            h_new[i][j] = omega * ((1 / a_p) * (
                    a_known_w * h_old[i][j - 1] + a_known_e * h_old[i][j + 1] + a_known_n * h_old[i - 1][j]
                    + a_known_s * h_old[i + 1][j] + source)) + ((1 - omega) * h_old[i][j])

    h_new[0][:] = h_new[1][:]
    h_new[:][n_x_max - 1] = h_new[:][n_x_max - 2]
    h_new[:][0] = h_new[:][1] + q * 2 * dx / 1000

    error = np.linalg.norm(h_new - h_old, 2)
    print('\nL2Norm = %0.5f' % error)

    print('iteration = ', iteration)

    plt.contourf(h_new)
    plt.gca().invert_yaxis()
    plt.axis('off')
    plt.grid()
    plt.colorbar().ax.set_ylabel('[m]')

    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # surf = ax.plot_surface(x, y, h_new)
    # fig.colorbar(surf)

    plt.pause(0.001)
    plt.show(block=False)
    plt.clf()

    # for i in range(n_y_max):
    #     for j in range(n_x_max):
    #         h_old[i][j] = h_new[i][j]

    h_old[:][:] = h_new[:][:]

print('L2Norm = ', np.linalg.norm(h_new))
print('iteration = ', iteration)

plt.contourf(h_new)
plt.gca().invert_yaxis()
plt.colorbar().ax.set_ylabel('[m]', rotation=270)
plt.savefig('Final_Result.png')
plt.show()

time2 = time.time()

print('\nTotal time = ', time2 - time1)

print()
