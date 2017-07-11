import random
import math
import matplotlib.pyplot as plt
import numpy
from scipy.interpolate import griddata
import sys
from scipy.stats import multivariate_normal



if __name__ == "__main__":
    #Generate gaussian image with noise N times
    if len(sys.argv) != 2:
        raise Exception("Syntax e.g.: python perlin.py 5")
    persistence = 0.9
    x, y = 256, 256 
    n = int(sys.argv[1])
    control = 0.01
    for i in range(n):
        noise = perlin2D(x,y,persistence)
        v1, v2 = numpy.mgrid[0:x, 0:y]
        v1 = numpy.array(v1, dtype=float)/127
        v2 = numpy.array(v2, dtype=float)/127
        pos = numpy.dstack((v1, v2))
        norm = multivariate_normal([0.5, 0.5], [[0.1, 0.0], [0.0, 0.1]])
        img = norm.pdf(pos)
        img = img/max(img.ravel())
        noise = noise/max(noise.ravel())
        signal = control*img
        noise = (1.0-control)*noise

        img = signal+noise
        print(numpy.average(img),numpy.std(img),)

        plt.imshow(signal+noise)
        plt.show()

    """
    #Generate a noise image:

    x, y = 800, 600 
    img = perlin2D(x,y)
    plt.imshow(img)
    plt.show()
    """
