import matplotlib
import matplotlib.pyplot as plt
import math
import random as rand

# creates a numerically simulated 2-D slice of a microstructure
# fine-grained assumes the rich phase form spheres separated by the eutectic phase
# directionally solidified assumes rectangular grains of equal size


def generate_canvas(sizeX, sizeY):
    # generates a 2-D canvas of dimensions sizeX pixels by sizeY pixels and fills it with 0's

    c = []
    for i in range(sizeX):
        c.append([])
        for j in range(sizeY):
            c[i].append(1)

    return c


def calibrate(canvas, length_fraction, r_um):
    # takes length fraction of sizeX of 1000um calibration bar to convert a length in microns to a length in pixels
    sizeX = len(canvas)
    cali_1000 = length_fraction * sizeX
    r_pix = r_um * cali_1000 / 1000
    return int(r_pix)


def generate_balls(canvas, volume_fraction, r_mean, r_std):
    # generates circles of rich phase and assigns them a value of 0 on the canvas
    # radius distribution is assumed to be gaussian
    # the circles may touch but may not overlap any amount

    sizeX = len(canvas)
    sizeY = len(canvas[0])
    canvas_volume = sizeX * sizeY
    v_ball = volume_fraction * canvas_volume
    n_balls = int(v_ball / (math.pi * r_mean**2))

    exclude = [[0, 0]]

    for i in range(n_balls):
        r = rand.gauss(r_mean, r_std)

        nuc_x = 0
        nuc_y = 0

        while [nuc_x, nuc_y] in exclude:
            nuc_x = rand.randint(0, sizeX)
            nuc_y = rand.randint(0, sizeY)

        else:
            for x in range(sizeX):
                for y in range(sizeY):
                    if (x-nuc_x)**2 + (y-nuc_y)**2 < r**2:
                        canvas[x][y] = 0

            for ex in range(sizeX):
                for ey in range(sizeY):
                    if (ex-nuc_x)**2 + (ey-nuc_y)**2 < (2*r)**2:
                        exclude.append([ex, ey])


def generate_file_fg(canvas, m):
    # outputs a file containing the canvas and returns the file
    sizeX = len(canvas)
    sizeY = len(canvas[0])

    file_out = str("/Volumes/Aroha Whanau/Bergman Research/USWS 2017/Microstructure Models/Sn-4pct-Zn FG " + str(sizeX) +
                   "x" + str(sizeY) + " x" + str(m) + ".txt")
    w = open(file_out, "w")

    for x in range(sizeX):
        for y in range(sizeY):
            w.write(str(canvas[x][y]))
            w.write("   ")
        w.write("\n")

    return file_out


def read_file(file):
    # reads a file containing the canvas and returns the canvas
    f = open(file)
    c = []
    for line in f:
        line_list = line.strip().split("   ")
        c.append(line_list)

    for x in range(len(c)):
        for y in range(len(c[0])):
            c[x][y] = int(c[x][y])

    return c


def generate_fg(x, y, m):
    # simulates a fine-grained microstructure
    canvas = generate_canvas(x, y)
    # canvas of our micrographs are 1392x1040

    if m == 20:
        # length fractions for calibrate function for various micrograph magnifications:
        # x20 magnification 0.0921
        # x40 magnification 0.206
        # x75 magnification 0.386
        r_mean = calibrate(canvas, 0.0921, 189)
        r_std = calibrate(canvas, 0.0921, 70)
    elif m == 40:
        r_mean = calibrate(canvas, 0.206, 189)
        r_std = calibrate(canvas, 0.206, 70)

    generate_balls(canvas, 0.667, r_mean, r_std)
    # volume fraction for Sn-4pct-Zn is 54.5%
    # area fraction, or 2/3 power of the volume fraction, is 66.7%

    file_in = generate_file_fg(canvas, m)

    canvas_w = read_file(file_in)

    plt.matshow(canvas_w, cmap="gray")
    plt.show()


# generate_fg(1000, 1000, 40)


def generate_file_ds(canvas, rh, m):
    # outputs a file containing the canvas and returns the file
    sizeX = len(canvas)
    sizeY = len(canvas[0])

    if rh == "r":
        file_out = str("/Volumes/Aroha Whanau/Bergman Research/USWS 2017/Microstructure Models/Sn-4pct-Zn DS " + str(sizeX) +
                       "x" + str(sizeY) + " x" + str(m) + " rectangles" + ".txt")
    elif rh == "h":
        file_out = str("/Volumes/Aroha Whanau/Bergman Research/USWS 2017/Microstructure Models/Sn-4pct-Zn DS " + str(sizeX) +
                       "x" + str(sizeY) + " x" + str(m) + " hexagons" + ".txt")

    w = open(file_out, "w")

    for x in range(sizeX):
        for y in range(sizeY):
            w.write(str(canvas[x][y]))
            w.write("   ")
        w.write("\n")

    return file_out


def generate_rectangles(canvas, volume_fraction, sizeP, sizeG):
    # generates rectangular grains with platelets that alternate rich and eutectic phases

    sizeX = len(canvas)
    sizeY = len(canvas[0])

    # volume fraction of the rich phase converted to a vertical offset of the cosine function
    vf_rich = -math.cos(volume_fraction * math.pi)

    # size of the platelets converted to a wavelength multiplier
    lamdbaP = sizeP / (2 * math.pi)

    nl_rectangles = int(sizeX / sizeG)

    grain_sizeX = int(sizeX / nl_rectangles)
    grain_sizeY = int(sizeY / nl_rectangles)

    for a in range(nl_rectangles):
        for b in range(nl_rectangles):
            m = rand.randint(0, 360)
            m_rad = m * math.pi / 180.0
            for x in range(grain_sizeX * a, grain_sizeX * (a + 1)):
                for y in range(grain_sizeY * b, grain_sizeY * (b + 1)):
                    if math.cos((x*math.cos(m_rad) - y*math.sin(m_rad)) / lamdbaP) + vf_rich > 0:
                        canvas[x][y] = 0
                    else:
                        canvas[x][y] = 1


def generate_hexagons(canvas, volume_fraction, sizeP, sizeG):
    # generates rectangular grains with platelets that alternate rich and eutectic phases

    sizeX = len(canvas)
    sizeY = len(canvas[0])

    # volume fraction of the rich phase converted to a vertical offset of the cosine function
    vf_rich = -math.cos(volume_fraction * math.pi)

    # size of the platelets converted to a wavelength multiplier
    lamdbaP = sizeP / (2 * math.pi)

    ul_hexagons = int(sizeX / sizeG)

    hex_height = int(sizeX / ul_hexagons)

    # hexagons in even spaces
    for a in range(ul_hexagons):
        for b in range(ul_hexagons):
            m = rand.randint(0, 360)
            m_rad = m * math.pi / 180.0
            for x in range(int(hex_height * a), int(hex_height * (a + 1/2))):
                for y in range(int((-1/math.sqrt(3)) * (x - hex_height * (a + 1/2)) + (2 / math.sqrt(3) * hex_height) * 1.5*b),
                               int(((1/math.sqrt(3)) * (x - hex_height * (a + 1/2))) + (2 / math.sqrt(3) * hex_height) * (1.5*b + 1))):
                    if not 0 < x < sizeX or not 0 < y < sizeY:
                        pass
                    elif math.cos((x*math.cos(m_rad) - y*math.sin(m_rad)) / lamdbaP) + vf_rich > 0:
                        canvas[x][y] = 0
                    else:
                        canvas[x][y] = 1
            for x in range(int(hex_height * (a + 1/2)), int(hex_height * (a + 1))):
                for y in range(int((1/math.sqrt(3)) * (x - hex_height * (a + 1/2)) + (2 / math.sqrt(3) * hex_height) * 1.5*b),
                               int(((-1/math.sqrt(3)) * (x - hex_height * (a + 1/2))) + (2 / math.sqrt(3) * hex_height) * (1.5*b + 1))):
                    if not 0 < x < sizeX or not 0 < y < sizeY:
                        pass
                    elif math.cos((x*math.cos(m_rad) - y*math.sin(m_rad)) / lamdbaP) + vf_rich > 0:
                        canvas[x][y] = 0
                    else:
                        canvas[x][y] = 1

    # hexagons in odd spaces
    for a in range(ul_hexagons + 1):
        for b in range(ul_hexagons):
            m = rand.randint(0, 360)
            m_rad = m * math.pi / 180.0
            for x in range(int(hex_height * (a - 1/2)), int(hex_height * a)):
                for y in range(int((-1 / math.sqrt(3)) * (x - hex_height * a) +
                                                           (2 / math.sqrt(3) * hex_height) * 1.5*b - 1.5 / math.sqrt(3) * hex_height),
                               int(((1 / math.sqrt(3)) * (x - hex_height * a)) +
                                                   (2 / math.sqrt(3) * hex_height) * (1.5*b + 1) - 1.5 / math.sqrt(3) * hex_height)):
                    if not 0 < x < sizeX or not 0 < y < sizeY:
                        pass
                    elif math.cos((x * math.cos(m_rad) - y * math.sin(m_rad)) / lamdbaP) + vf_rich > 0:
                        canvas[x][y] = 0
                    else:
                        canvas[x][y] = 1
            for x in range(int(hex_height * a), int(hex_height * (a + 1/2))):
                for y in range(int((1 / math.sqrt(3)) * (x - hex_height * a) +
                                           ((2 / math.sqrt(3) * hex_height) * 1.5*b) - 1.5 / math.sqrt(3) * hex_height),
                               int(((-1 / math.sqrt(3)) * (x - hex_height * a)) +
                                            (2 / math.sqrt(3) * hex_height) * (1.5*b + 1) - 1.5 / math.sqrt(3) * hex_height)):
                    if not 0 < x < sizeX or not 0 < y < sizeY:
                        pass
                    elif math.cos((x * math.cos(m_rad) - y * math.sin(m_rad)) / lamdbaP) + vf_rich > 0:
                        canvas[x][y] = 0
                    else:
                        canvas[x][y] = 1


def generate_ds(x, y, m):
    # simulates a directionally-solidified microstructure
    canvas = generate_canvas(x, y)

    if m == 20:
        sizeP = calibrate(canvas, 0.0921, 300)
        sizeG = calibrate(canvas, 0.0921, 1500)
    elif m == 40:
        sizeP = calibrate(canvas, 0.206, 300)
        sizeG = calibrate(canvas, 0.206, 1500)

    generate_rectangles(canvas, 0.9, sizeP, sizeG)
    file_in_r = generate_file_ds(canvas, "r", m)
    #
    # generate_hexagons(canvas, 0.545, sizeP, sizeG)
    # file_in_h = generate_file_ds(canvas, "h", m)

    canvas_w = read_file(file_in_r)

    plt.matshow(canvas_w, cmap="gray")
    plt.show()


generate_ds(500, 500, 40)
