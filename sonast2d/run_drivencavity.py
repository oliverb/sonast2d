import cProfile
import numpy
from math import pi
from nast2d import solver

class DrivenCavity:
    def __init__(self, imax=64, jmax=64):
        self.xlength = 1.0
        self.ylength = 1.0
        self.t_end = 20.0

        self.imax = imax
        self.jmax = jmax 

        self.Re = 50.0
        self.omega = 2.0-2*pi*self.xlength/self.imax
        self.tau = 0.5
        self.alpha = 0.9
        self.eps = 1E-7
        self.max_it = 2000

        self.boundary = {'scenario':'cavity'}

        self.geo = numpy.zeros((imax+2,jmax+2), dtype=numpy.int32)

# cProfile.run('solver.solve(solver.DrivenCavity())')
solver.solve(DrivenCavity(128,128))
