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

# cvs_res = int(input('Canvas resolution (in pixels)? '))
# cvs_mag = int(input('Choose magnification: x20, x40, or x75? '))
cvs_res = 1000
cvs_mag = 40

# Initialising values
x = np.arange(0, cvs_res, 1)
y = np.arange(0, cvs_res, 1)
Xg, Yg = np.meshgrid(x, y, sparse=False)
Vg = np.zeros((cvs_res, cvs_res))
Vtest = np.zeros((cvs_res, cvs_res))
r = 0
m = 0
a_prev = 0
a_tot = 0

r_mean_px = calibrate(r_mean, cvs_mag, cvs_res)
r_std_px = calibrate(r_std, cvs_mag, cvs_res)
overlap = 0.7               # defines percentage of the grain area not allowed to overlap with neighbouring grains

# defines allowed nucleation sites for grains, randomly picks the center of each circular grain from these
site_x = np.arange(r_mean_px//2, cvs_res-r_mean_px//2, int(r_mean_px))
site_y = np.arange(r_mean_px//2, cvs_res-r_mean_px//2, int(r_mean_px))
lng = len(site_x) - 1
Xs, Ys = np.meshgrid(site_x, site_y, sparse=False)
np.random.shuffle(Xs.flat)
np.random.shuffle(Ys.flat)

while a_tot < area_frac * cvs_res**2:
    while a_tot <= a_prev + overlap * (np.pi * r_mean_px**2):
        r = rand.gauss(r_mean_px, r_std_px)
        c = 2 * r / 2.35482

        # allows fluctuations around selected nucleation site to recover some randomness
        nuc_x = Xs[m%lng, m//lng] + (rand.randint(0, 2*r_mean_px)-r_mean_px)
        nuc_y = Ys[m%lng, m//lng] + (rand.randint(0, 2*r_mean_px)-r_mean_px)
        print('Testing grain #' + str(m+1) + ' at point: (' + str(nuc_x) + ', ' + str(nuc_y) + ')')

        V = np.exp(-(((Xg-nuc_x)**2)/(2*(c**2)) + ((Yg-nuc_y)**2)/(2*(c**2))))
        V_thresh = np.where(V >= 0.5, 1, 0)

        Vtest = Vg + V_thresh
        Ztest = np.where(Vtest >= 0.5, 1, 0)
        a_tot = np.sum(Ztest)

    m += 1
    print('Drawing grain #' + str(m) + ' at point: (' + str(nuc_x) + ', ' + str(nuc_y) + ')')

    Vg += V_thresh
    Vtest = Vg
    a_prev = a_tot

Z = np.where(Vg >= 0.5, 1, 0)
a_frac = np.sum(Z) / cvs_res**2
print('Area fraction:', a_frac)

end = time.time()
comp_time = end - start
print('Total computation time (in seconds):', comp_time)

plt.matshow(Z, cmap="gray")
plt.show()
