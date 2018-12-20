"""
Doc string goes here ~ ~ ~
"""

import time
import math
import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
import random as rand
from Constants import *

start = time.time()

cvs_res = 500

x = np.arange(0, cvs_res, 1)
y = np.arange(0, cvs_res, 1)
Xg, Yg = np.meshgrid(x, y, sparse=True)
Vg = np.zeros((cvs_res, cvs_res))

r_mean_px = calibrate(r_mean, 40, 500)
r_std_px = calibrate(r_std, 40, 500)
n = int((area_frac * cvs_res**2) / (math.pi * r_mean_px**2))
print('Total number of grains:', n)
m = 0

exc = [[0, 0]]
for i in range(n):
    r = rand.gauss(r_mean_px, r_std_px)
    c = 2*r/2.35482

    while True:
        nuc_x = rand.randint(0, cvs_res)
        nuc_y = rand.randint(0, cvs_res)
        min_dist = cvs_res * 2
        for pt in exc:
            dist = linalg.norm(np.array([nuc_x, nuc_y]) - np.array(pt))
            if dist < min_dist:
                min_dist = dist

        if min_dist >= 2*r_mean_px:
            m += 1
            print('Drawing grain #' + str(m) + ' at point: (' + str(nuc_x) + ', ' + str(nuc_y) + ')')
            exc.append([nuc_x, nuc_y])
            break

    Vg = Vg + np.exp(-((Xg-nuc_x)**2 + (Yg-nuc_y)**2)/(2*c)**2)

Z = np.where(Vg >= 0.5, 1, 0)

end = time.time()
comp_time = start - end
print('Total computation time:', comp_time)

plt.matshow(Z, cmap="gray")
plt.show()
