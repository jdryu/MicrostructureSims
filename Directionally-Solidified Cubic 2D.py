"""
Doc string goes here ~ ~ ~
"""

import time
import numpy as np
import scipy.spatial as sp
import matplotlib.pyplot as plt
from Constants import *

start = time.time()

# cvs_res = int(input('Canvas resolution (in pixels)? '))
# cvs_mag = int(input('Choose magnification: x20, x40, or x75? '))
cvs_res = 1000
cvs_mag = 40
name = user_dir + date + ' DSC 2D x' + str(cvs_mag)

# Initialising values
x = np.arange(0, cvs_res, 1)
y = np.arange(0, cvs_res, 1)
Xg, Yg = np.meshgrid(x, y, sparse=False)
Vg = np.zeros((cvs_res, cvs_res))

sizeP_px = calibrate(sizeP, cvs_mag, cvs_res)
sizeG_px = calibrate(sizeG, cvs_mag, cvs_res)
cx = 2 * 0.4*sizeP_px / 2.35482
cy = 2 * 0.1*sizeP_px / 2.35482

# defines allowed nucleation sites for grains
Gsite_x = np.arange(sizeG_px//2, cvs_res-sizeG_px//2, sizeG_px)
Gsite_y = np.arange(sizeG_px//2, cvs_res-sizeG_px//2, sizeG_px)
Xs, Ys = np.meshgrid(Gsite_x, Gsite_y, sparse=False)
Gpoints = np.column_stack((Xs.flat, Ys.flat))
l, w = Gpoints.shape
perturb = np.random.randn(l, w) * sizeG_px//4               # allows fluctuations around nucleation site to recover some randomness
Gpoints = Gpoints + perturb

# define nucleation sites for microstructure
Psite_x = np.arange(sizeP_px//2, cvs_res, sizeP_px)
Psite_y = np.arange(sizeP_px//2, cvs_res, sizeP_px)
Xp, Yp = np.meshgrid(Psite_x, Psite_y, sparse=False)
Ppoints = np.column_stack((Xp.flat, Yp.flat))

for pt in Ppoints:
    midx = sp.distance.cdist([pt], Gpoints).argmin()
    # plt.plot(pt[0], pt[1], 'o', markersize=2, color=colorlist[midx%len(colorlist)])
    Vx = np.exp(-(((Xg - pt[0]) ** 2) / (2 * (cx ** 2)) + ((Yg - pt[1]) ** 2) / (2 * (cy ** 2))))
    Vy = np.exp(-(((Xg - pt[0]) ** 2) / (2 * (cy ** 2)) + ((Yg - pt[1]) ** 2) / (2 * (cx ** 2))))
    V_thresh = np.where(Vx+Vy >= 0.5, 1, 0)
    Vg += V_thresh


# plt.plot(Gpoints[:, 0], Gpoints[:, 1], 'ko', markersize=5)
# plt.axis([0, cvs_res, 0, cvs_res])
# plt.show()


Z = np.where(Vg >= 0.5, 1, 0)
a_frac = np.sum(Z) / cvs_res**2
print('Area fraction:', a_frac)

end = time.time()
comp_time = end - start
print('Total computation time (in seconds):', comp_time)

plt.matshow(Z, cmap='gray')
plt.savefig(name + '.png')
plt.show()

np.save(name + '.npy', Z)
