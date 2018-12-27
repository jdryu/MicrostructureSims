"""
Doc string goes here ~ ~ ~
"""

import os
import datetime


def calibrate(r_um, mag, res):
    # takes relative length of 1000um calibration bar to convert a length in microns to a length in pixels
    if mag == 20:
        len_frac = 0.0921
    elif mag == 40:
        len_frac = 0.206
    elif mag == 75:
        len_frac = 0.386
    else:
        # default to x20 magnification
        len_frac = 0.0921
    cali_1000 = len_frac * res
    r_pix = r_um * cali_1000 / 1000
    return int(r_pix)


date = str(datetime.date.today())
user_dir = os.getcwd() + "/Simulated Models/"

r_mean = 189        # mean radius of equiaxed grains (approx as spheres)
r_std = 70          # standard deviation of the radius of equiaxed grains

sizeP = 300         # approximate mean spacing between Zn-Sn platelets
sizeG = 1500        # approximate mean size of directionally solidified grains

vol_frac = 0.545        # pure-to-eutectic volume fraction for Sn-4pct-Zn
area_frac = 0.667       # pure-to-eutectic area fraction for Sn-4pct-Zn
