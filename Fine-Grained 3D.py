"""
Doc string goes here ~ ~ ~
"""

import time
import numpy as np
import matplotlib.pyplot as plt
import random as rand
from Constants import *

start = time.time()

# cvs_res = int(input('Canvas resolution (in pixels)? '))
# cvs_mag = int(input('Choose magnification: x20, x40, or x75? '))
cvs_res = 500
cvs_mag = 75
name = user_dir + date + ' FG 3D x' + str(cvs_mag)

# Initialising values
x = np.arange(0, cvs_res, 1)
y = np.arange(0, cvs_res, 1)
z = np.arange(0, cvs_res, 1)
Xg, Yg, Zg = np.meshgrid(x, y, z, sparse=False)
Vg = np.zeros((cvs_res, cvs_res, cvs_res))
Vtest = np.zeros((cvs_res, cvs_res, cvs_res))
tcount = 0
v_prev = 0
v_tot = 0
m = 0

r_mean_px = calibrate(r_mean, cvs_mag, cvs_res)
r_std_px = calibrate(r_std, cvs_mag, cvs_res)
overlap = 0.9               # defines percentage of the sphere volume not allowed to overlap with neighbouring sphere
print('Approximate number of grains:', vol_frac*cvs_res**3//(4/3*np.pi * r_mean_px**3))

# defines allowed nucleation sites for grains, randomly picks the center of each spherical grain from these
site_x = np.arange(r_mean_px//2, cvs_res-r_mean_px//2, int(r_mean_px))
site_y = np.arange(r_mean_px//2, cvs_res-r_mean_px//2, int(r_mean_px))
site_z = np.arange(r_mean_px//2, cvs_res-r_mean_px//2, int(r_mean_px))
lng = len(site_x) - 1
Xs, Ys, Zs = np.meshgrid(site_x, site_y, site_z, sparse=False)
np.random.shuffle(Xs.flat)
np.random.shuffle(Ys.flat)
np.random.shuffle(Zs.flat)

while v_tot < vol_frac * cvs_res**3:
    while True:
        r = rand.gauss(r_mean_px, r_std_px)
        c = 2 * r / 2.35482

        # allows fluctuations around selected nucleation site to recover some randomness
        nuc_x = Xs[m%lng, (m//lng)%lng, (m//lng**2)%lng] + (rand.randint(0, 2*r_mean_px)-r_mean_px)
        nuc_y = Ys[m%lng, (m//lng)%lng, (m//lng**2)%lng] + (rand.randint(0, 2*r_mean_px)-r_mean_px)
        nuc_z = Zs[m%lng, (m//lng)%lng, (m//lng**2)%lng] + (rand.randint(0, 2*r_mean_px)-r_mean_px)
        print('Testing grain #' + str(m+1) + ' at point: (' + str(nuc_x) + ', ' + str(nuc_y) + ', ' + str(nuc_z) + ')')

        V = np.exp(-(((Xg-nuc_x)**2)/(2*(c**2)) + ((Yg-nuc_y)**2)/(2*(c**2)) + ((Zg-nuc_z)**2)/(2*(c**2))))
        V_thresh = np.where(V >= 0.5, 1, 0)

        Vtest = Vg + V_thresh
        Wtest = np.where(Vtest >= 0.5, 1, 0)
        v_tot = np.sum(Wtest)

        if v_tot >= v_prev + overlap * (4/3 * np.pi * r ** 3):
            tcount = 0
            break
        elif tcount > 10:
            print('Failed too many tries to fit without overlap. Using latest attempt for nucleation site.')
            tcount = 0
            break
        else:
            tcount += 1

    print('Drawing grain #' + str(m+1) + ' at point: (' + str(nuc_x) + ', ' + str(nuc_y) + ', ' + str(nuc_z) + ')')
    print('Current volume fraction at: ', v_tot/cvs_res**3)
    print('Current computation time (in seconds) at: ', time.time()-start)
    m += 1
    Vg += V_thresh
    Vtest = Vg
    v_prev = v_tot

W = np.where(Vg >= 0.5, 1, 0)
v_frac = np.sum(W) / cvs_res**3
print('Volume fraction:', v_frac)

end = time.time()
comp_time = end - start
print('Total computation time (in seconds):', comp_time)

plt.matshow(W, cmap="gray")
plt.savefig(name + '.png')
plt.show()

np.save(name + '.npy', W)
