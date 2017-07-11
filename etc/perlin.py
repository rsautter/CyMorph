import random
import math
import matplotlib.pyplot as plt
import numpy
from scipy.interpolate import griddata
import sys
from scipy.stats import multivariate_normal,signaltonoise


def perlin2D(x, y, persistence = 0.7):
    """
	Adapated from: http://code.activestate.com/recipes/578470-perlin-noise-generator/
	
    """
    imgAr = [[0.0 for i in range(x)] for j in range(y)]
    totAmp = 0.0
    octaves = int(math.log(max(x, y), 2.0))
    for k in range(octaves):
        freq = 2 ** k
        amp = persistence ** k
        totAmp += amp
        n = freq + 1; m = freq + 1 
        ar = [[random.random() * amp for i in range(n)] for j in range(m)]
        nx = x / (n - 1.0); ny = y / (m - 1.0)

        for ky in range(y):
            for kx in range(x):
                i = int(kx / nx)
                j = int(ky / ny)
                
                # distances
                dx0 = kx - i * nx
                dx1 = nx - dx0
                dy0 = ky - j * ny
                dy1 = ny - dy0

                #interpolating (bilinear)
                z = ar[j][i] * dx1 * dy1
                z += ar[j][i + 1] * dx0 * dy1
                z += ar[j + 1][i] * dx1 * dy0
                z += ar[j + 1][i + 1] * dx0 * dy0
                z /= nx * ny 
                imgAr[ky][kx] += z 
    return numpy.array(imgAr)

if __name__ == "__main__":
    #Generate gaussian image with noise N times
    if len(sys.argv) != 2:
        raise Exception("Syntax e.g.: python perlin.py 5")
    persistence = 0.9
    x, y = 256, 256 
    n = int(sys.argv[1])
    control = 0.75
    for i in range(n):
        noise = perlin2D(x,y,persistence)
        v1, v2 = numpy.mgrid[0:x, 0:y]
        v1 = numpy.array(v1, dtype=float)/float(x)
        v2 = numpy.array(v2, dtype=float)/float(y)
        pos = numpy.dstack((v1, v2))
        norm = multivariate_normal([0.5, 0.5], [[0.1, 0.0], [0.0, 0.1]])
        img = norm.pdf(pos)
        img = img/max(img.ravel())
        noise = noise/max(noise.ravel())
        signal = control*img
        noise = (1.0-control)*noise

        img = signal+noise
        print(signaltonoise(img,axis = None)[0])

        plt.imshow(signal+noise)
        plt.show()

    """
    #Generate a noise image:

    x, y = 800, 600 
    img = perlin2D(x,y)
    plt.imshow(img)
    plt.show()
    """
