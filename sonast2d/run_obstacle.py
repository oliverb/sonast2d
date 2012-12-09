from math import pi

import numpy

from nast2d import solver
from splines import bspline

def gen_geo(imax, jmax, xlength ,ylength, indicator):
    dx = xlength/imax
    dy = ylength/jmax

    geo = numpy.zeros((imax+2,jmax+2), dtype=numpy.int32)

    for j in range(1, jmax+1):
        for i in range(1, imax+1):
            geo[i,j] = indicator(i*dx, j*dy)

    return geo

def obstacle_chi(c):
    curve = bspline.BasisSpline(3, 
                numpy.array([0.0, 4.0/3.0, 8.0/3.0, 4.0]))

    def chi(x, y):
        if x <= 2.0:
            return 0
        if x >= 6.0:
            return 0

        sval = curve.evaluate(c, x-2.0)

        if abs(y - 4.0) <= sval:
            return 1
        else:
            return 0

    return chi



class ObstacleFlow:
    def __init__(self, c, imax=192, jmax=64):
        self.xlength = 24.0
        self.ylength = 8.0
        self.Re = 80.0
        self.t_end = 100.0

        self.boundary = {'scenario':'obstacle',
                         'north':'slip',
                         'south':'slip',
                         'east':'slip',
                         'west':'slip'
                        }

        self.imax = imax
        self.jmax = jmax

        h = self.ylength/self.jmax
        self.omega = 2.0-2*pi*h
        self.tau = 0.5
        self.alpha = 0.9
        self.eps = 1E-7
        self.max_it = 5000


        # Generate geo
        indicator = obstacle_chi(c)
        self.geo = gen_geo(self.imax, self.jmax, self.xlength,
                           self.ylength, indicator)

c = numpy.array([0.0, 0.8, 0.979795, 0.979795, 0.8, 0.0])
problem = ObstacleFlow(c)

solver.solve(problem)