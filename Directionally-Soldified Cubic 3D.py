"""
Doc string goes here ~ ~ ~
"""

import time
import numpy as np
import scipy.spatial as sp
from mayavi import mlab
from Constants import *

start = time.time()

# cvs_res = int(input('Canvas resolution (in pixels)? '))
# cvs_mag = int(input('Choose magnification: x20, x40, or x75? '))
cvs_res = 500
cvs_mag = 40
name = user_dir + date + ' DSC 3D x' + str(cvs_mag)

# Initialising values
x = np.arange(0, cvs_res, 1)
y = np.arange(0, cvs_res, 1)
z = np.arange(0, cvs_res, 1)
Xg, Yg = np.meshgrid(x, y, sparse=False)
Vg = np.zeros((cvs_res, cvs_res))

sizeP_px = calibrate(sizeP, cvs_mag, cvs_res)
sizeG_px = calibrate(sizeG, cvs_mag, cvs_res)
sx2 = (2 * 0.7*sizeP_px / 2.35482)**2
sy2 = (2 * 0.1*sizeP_px / 2.35482)**2

# defines allowed nucleation sites for grains
Gsite_x = np.arange(sizeG_px//2, cvs_res-sizeG_px//2, sizeG_px)
Gsite_y = np.arange(sizeG_px//2, cvs_res-sizeG_px//2, sizeG_px)
Xs, Ys = np.meshgrid(Gsite_x, Gsite_y, sparse=False)
Gpoints = np.column_stack((Xs.flat, Ys.flat))
l, w = Gpoints.shape
perturb = np.random.randn(l, w) * sizeG_px//4               # allows fluctuations around nucleation site to recover some randomness
Gpoints = Gpoints + perturb

# define properties for each grain
theta = np.random.rand(l) * 2 * np.pi
phi = np.random.rand(l) * np.pi/36 - np.pi/18
ax = np.cos(theta)**2 / (2 * sx2) + np.sin(theta)**2 / (2 * sy2)
bx = -np.sin(2*theta) / (4 * sx2) + np.sin(2*theta) / (4 * sy2)
cx = np.sin(theta)**2 / (2 * sx2) + np.cos(theta)**2 / (2 * sy2)
ay = np.cos(theta)**2 / (2 * sy2) + np.sin(theta) ** 2 / (2 * sx2)
by = -np.sin(2*theta) / (4 * sy2) + np.sin(2*theta) / (4 * sx2)
cy = np.sin(theta)**2 / (2 * sy2) + np.cos(theta) ** 2 / (2 * sx2)

# define nucleation sites for microstructure
Psite_x = np.arange(-sizeP_px//2, cvs_res+sizeP_px//2, sizeP_px)
Psite_y = np.arange(-sizeP_px//2, cvs_res+sizeP_px//2, sizeP_px)
Xp, Yp = np.meshgrid(Psite_x, Psite_y, sparse=False)
Ppoints = np.column_stack((Xp.flat, Yp.flat))

for pt in Ppoints:
    midx = sp.distance.cdist([pt], Gpoints).argmin()
    ptx = (np.cos(phi[midx]) * (pt[0]-Gpoints[midx][0]) - np.sin(phi[midx]) * (pt[1]-Gpoints[midx][1])) + Gpoints[midx][0]
    pty = (np.sin(phi[midx]) * (pt[0]-Gpoints[midx][0]) + np.cos(phi[midx]) * (pt[1]-Gpoints[midx][1])) + Gpoints[midx][1]

    px = ax[midx]*((Xg-ptx)**2) + 2*bx[midx]*(Xg-ptx)*(Yg-pty) + cx[midx]*((Yg-pty)**2)
    py = ay[midx]*((Xg-ptx)**2) + 2*by[midx]*(Xg-ptx)*(Yg-pty) + cy[midx]*((Yg-pty)**2)

    Vx = np.exp(-px)
    Vy = np.exp(-py)
    Vg += Vx + Vy

Z = np.where(Vg >= 0.5, 1, 0)
W = np.repeat(Z[:, :, np.newaxis], cvs_res, axis=2)
np.save(name + '.npy', W)
v_frac = np.sum(W) / cvs_res**3
print('Volume fraction:', v_frac)

mlab.figure()
mlab.contour3d(W)

end = time.time()
comp_time = end - start
print('Total computation time (in seconds):', comp_time)

mlab.show()
mlab.savefig(name + '.png')
