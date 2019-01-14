import numpy as np
import matplotlib.pyplot as plt


def find_neighbors(xc, yc, grid_size):
    xround = int(round(xc))
    yround = int(round(yc))

    up_left = ((xround - 1) % grid_size, yround % grid_size)
    up_right = (xround % grid_size, yround % grid_size)
    low_left = ((xround - 1) % grid_size, (yround - 1) % grid_size)
    low_right = (xround % grid_size, (yround - 1) % grid_size)
    return up_left, up_right, low_left, low_right


def calc_upper_left_area(xc, yc):
    return (round(xc) - (xc - 0.5)) * ((yc + 0.5) - round(yc))


def calc_upper_right_area(xc, yc):
    return ((xc + 0.5) - round(xc)) * ((yc + 0.5) - round(yc))


def calc_lower_left_area(xc, yc):
    return (round(xc) - (xc - 0.5)) * (round(yc) - (yc - 0.5))


def calc_lower_right_area(xc, yc):
    return ((xc + 0.5) - round(xc)) * (round(yc) - (yc - 0.5))


def area_weight(xc, yc, grid, grid_size):
    up_left, up_right, low_left, low_right = find_neighbors(xc, yc, grid_size)
    grid[up_left] += calc_upper_left_area(xc, yc)
    grid[up_right] += calc_upper_right_area(xc, yc)
    grid[low_left] += calc_lower_left_area(xc, yc)
    grid[low_right] += calc_lower_right_area(xc, yc)
    return


n = 10  # max range of data points
num_points = 50  # only for random number of data points to appear
weighted_grid = np.zeros((n, n))
rand_data_x = np.random.rand(num_points) * (n - 1)
rand_data_y = np.random.rand(num_points) * (n - 1)
# rand_data_x = np.ones(num_points) * 0
# rand_data_y = np.ones(num_points) * 0

for i in xrange(num_points):
    xc = rand_data_x[i]
    yc = rand_data_y[i]
    area_weight(xc, yc, weighted_grid, n)

plt.figure(1)
plt.contourf(weighted_grid)
# plt.show()
# plt.figure(2)
plt.scatter(rand_data_x, rand_data_y, c='r')
# plt.colorbar()
plt.show()
# plt.figure(3)
# plt.contourf(rand_data_x, rand_data_y)
# plt.show()
