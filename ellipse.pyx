from libc.math cimport pow, sqrt, cos, sin
from math import radians

import galaxyIO
from libc.math cimport pow, sqrt, cos, sin
from math import radians

cimport cython 

cdef class ellipse:
    cdef public float rp, id
    cdef public float angle, petrosianMag, sky, maxRad, minRad, fwhm
    cdef public float posx, posy

    # Warning! - this builds an empty object
    cdef void setEmpty(self):
        print("Warning! - Empty Ellipse Object")
        self.rp = 0.0
        self.id = 0.0
        self.posx = 0.0
        self.posy = 0.0
        self.fwhm = 0.0
        self.petrosianMag = 0.0
        self.sky = 0.0
        self.angle = 0.0
        self.maxRad = 0.0
        self.minRad = 0.0

    def __cinit__(self, dic=None, data=None,calibratedData=None,float sky=-1.0):
        cdef float angulo,scaleA,scaleB
        if(dic is None or data is None or calibratedData is None):
            self.setEmpty()
        else:
            self.rp = float(galaxyIO.getFirstValue(dic, data, 'PETRO_RADIUS'))
            self.id = float(galaxyIO.getFirstValue(dic, data, 'NUMBER'))
            self.posx = float(galaxyIO.getFirstValue(dic, data, 'X_IMAGE'))
            self.posy = float(galaxyIO.getFirstValue(dic, data, 'Y_IMAGE'))
            angulo = radians(float(galaxyIO.getFirstValue(dic, data, 'THETA_IMAGE')))
            scaleA = float(galaxyIO.getFirstValue(dic, data, 'A_IMAGE'))
            scaleB = float(galaxyIO.getFirstValue(dic, data, 'B_IMAGE'))
            self.petrosianMag = float(galaxyIO.getFirstValue(dic, data, 'FLUX_PETRO'))
            self.fwhm = float(galaxyIO.getFirstValue(dic, data, 'FWHM_IMAGE'))
            if(sky == -1.0):
                self.sky = float(galaxyIO.getFirstValue(dic, data, 'BACKGROUND'))
            else:
                self.sky = sky 
            if (calibratedData):
                self.petrosianMag = pow(10.0,-0.4*self.petrosianMag)
            #self.pos = (x,y)
            self.angle = angulo
            self.maxRad = scaleA*self.rp
            self.minRad = scaleB*self.rp
 
    ##Encontra a escala da elipse que corta x, y
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef float findScale(self, float x,float  y):
        cdef float dx = x-float(self.posx)
        cdef float dy = y-float(self.posy)
        return sqrt(pow((dx*cos(self.angle)+dy*sin(self.angle))/(self.maxRad),2.0) + \
                pow((dx*sin(self.angle)-dy*cos(self.angle))/(self.minRad),2.0))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef float getCartesianDist(self, float x, float y):
        cdef float dx = x-float(self.posx)
        cdef float dy = y-float(self.posy)
        return (dx**2.0+dy**2.0)**0.5

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef tuple findPoint(self, float scale,float rho):
        x, y = scale*self.maxRad*cos(rho), scale*self.minRad*sin(rho)
        x2, y2 = self.posx+x*cos(self.angle)-y*sin(self.angle), self.posy+x*sin(self.angle)+y*cos(self.angle)
        return x2, y2
