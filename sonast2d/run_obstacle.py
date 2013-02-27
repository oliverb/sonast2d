#!/usr/bin/env python
from math import pi
import sys
import os

import numpy
import scipy.optimize

from nast2d import solver
from splines import bspline
from tools import numpy2hopspack

def gen_geo(imax, jmax, xlength ,ylength, indicator):
    dx = xlength/imax
    dy = ylength/jmax

    geo = numpy.zeros((imax+2,jmax+2), dtype=numpy.int32)

    for j in range(1, jmax+1):
        for i in range(1, imax+1):
            geo[i,j] = indicator((i-1.)*dx+.5*dx, (j-1.)*dy+.5*dy)

    return geo

# c = coefficients
def obstacle_chi(c):
    N, = c.shape

    # For N degrees of freedom we need N-2 knots
    curve = bspline.BasisSpline(3, numpy.linspace(0.0, 4.0, N-2))

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
    def __init__(self, c, imax=3*32, jmax=32, name="obstacle"):
        self.xlength = 24.0
        self.ylength = 8.0
        self.Re = 10.0
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
        self.omega = 1.7 #2.0-2*pi*h
        self.tau = 0.5
        self.alpha = 0.9
        self.eps = 1E-10
        self.max_it = 5000

        self.name = name


        # Generate geo
        indicator = obstacle_chi(c)
        self.geo = gen_geo(self.imax, self.jmax, self.xlength,
                           self.ylength, indicator)

# Routines to calculate nice circulat starting values ...
def half_circle(x):
    if x >= 2.0:
        return 0.0
    x = x-1.0
    return numpy.sqrt(1.1-x*x)

# N -> Number of unknowns
def nice_coeff(N, border='yes'):
    curve = bspline.BasisSpline(3, numpy.linspace(0.0, 4.0, N-2))

    xvals = numpy.linspace(0.0, 4.0, N)

    # Matrix with values of basis functions at xvals
    # so that B*c gives a vector with the evaluations of the spline with
    # coefficients c at each x \in xvals
    B = curve.compute_basis_matrix(xvals)

    # Prepare vector with desired values at xvals
    f = []
    for x in xvals:
        f.append(half_circle(x))
    
    f = numpy.array(f)

    if border=='no':
        f = f[1:N-1]
        B = B[1:N-1,1:N-1]
        Bineq = curve.compute_basis_matrix(numpy.linspace(0.0, 4.0, N-1))
        Bineq = Bineq[1:N-2,1:N-1]

    def grad_func(c):
        return 0.5*numpy.dot(c, numpy.dot(B, c)) - numpy.dot(c, f)

    def dgrad_func(c):
        return numpy.dot(B, c) - f

    def ineq_func(c):
        if border=='no':
            return numpy.dot(Bineq, c)
        numpy.dot(B, c)

    mini = scipy.optimize.fmin_slsqp(grad_func, f)
    
    return mini

def apps_routine(imax, jmax):
    ifilename = sys.argv[1]
    ofilename = sys.argv[2]
    tag = sys.argv[3]

    Xpart = numpy.genfromtxt(ifilename)[1:]

    X = numpy.r_[0.0, Xpart, 0.0]

    problem = ObstacleFlow(X, name="foo"+tag, imax=imax, jmax=jmax)
    max_diss = solver.solve(problem)

    with open(ofilename, "w") as ofile:
        print >>ofile, "%.12f"%max_diss

def hops_config():
    N = int(sys.argv[3])
    x0 = nice_coeff(N, border='no')
    curve= bspline.BasisSpline(3, numpy.linspace(0.0, 4.0, N-2))

    # Ineq. constraints
    lowvec = numpy.r_[numpy.zeros(N-1), pi/2.0]
    upvec = numpy.r_[numpy.ones(N-1)*2.0, pi]

    xvals = numpy.linspace(0.0, 4.0, N-1)
    B = curve.compute_basis_matrix(xvals)
    I = curve.compute_basis_integrals()
    A = numpy.vstack([B, I])

    configfile = sys.argv[2]

    cur_file = os.path.abspath(__file__)
    numpy2hopspack.makecfg(N, cur_file, configfile, x0=x0, Aineq=A, lineq=lowvec,
                            uineq=upvec, P=8, gss_step_top=0.01)
    
def apps_config():
    N = int(sys.argv[3])

    x0 = nice_coeff(N, border='no')
    curve = bspline.BasisSpline(3, numpy.linspace(0.0, 4.0, N-2))

    with open(sys.argv[2], "w") as cf:
        cf.write('# SAMPLE APPSPACK INPUT FILE\n')
        cf.write('@ "Linear"\n')
#        cf.write('"Upper" vector 4  10 10 10 10\n')
#        cf.write('"Lower" vector 4  -10 0 0 0\n')
#        cf.write('"Scaling" vector 4 1 1 1 1\n')
#        cf.write('"Inequality Upper" vector 2    DNE   -1\n')
#        cf.write('"Inequality Lower" vector 2    -10  DNE\n')
#        cf.write('"Equality Bound" vector 1 3\n')
#        cf.write('"Inequality Matrix"  matrix  2 4 \n')
#        cf.write('-1 -1 -1 -1\n')
#        cf.write(' 1 -1  1 -1\n')
#        cf.write('"Equality Matrix" matrix 1 4\n')
#        cf.write('2 0 2 -7\n')

        # Upper

        # Scaling vector
        cf.write('"Scaling" vector %d ' % (N-2))
        for i in range(1, N-1):
            cf.write('1.0 ')
        cf.write('\n')

        # Pointwise inequalities + Volume bounds
        cf.write('"Inequality Upper" vector %d ' % (N-2))
        for i in range(1, N-2):
            cf.write('2.0 ')
        cf.write('%f ' % pi)
        cf.write('\n')

        cf.write('"Inequality Lower" vector %d ' % (N-2))
        for i in range(1, N-2):
            cf.write('0.0 ')
        cf.write('%f ' % (pi/2.0))
        cf.write('\n')        

        xvals = numpy.linspace(0.0, 4.0, N-1)
        B = curve.compute_basis_matrix(xvals)
        I = curve.compute_basis_integrals()

        cf.write('"Inequality Matrix" matrix %d %d \n' % (N-2, N-2))
        for i in range(1, N-2):
            for j in range(1, N-1):
                cf.write('%f ' % B[i,j])
            cf.write('\n')
        for j in range(1, N-1):
            cf.write('%f ' % I[j])
        cf.write('\n')



        cf.write('@@\n')
        cf.write('@ "Evaluator"\n')
        cf.write('"Executable Name" string "run_obstacle.py"\n')
        cf.write('"Input Prefix" string "obstacle_input"\n')
        cf.write('"Output Prefix" string "obstacle_output"\n')
        cf.write('@@\n')
        cf.write('@ "Solver" \n')
        cf.write('"Debug" int 3\n')

        # Initial vector
        cf.write('"Initial X" vector %d ' % (N-2))
        for i in range(0, N-2):
            cf.write("%f " % x0[i])
        cf.write('\n')

        cf.write('"Step Tolerance" double 1e-3\n')
        cf.write('@@\n')
        print numpy.dot(I[1:N-1], x0)
        #print numpy.dot(B[1:N-1,1:N-1], x0)
        print str(pi/2.0)


def main():
    # c = numpy.array([0.0, 0.8, 0.979795, 0.979795, 0.8, 0.0])
    # c = nice_coeff(8)
    # problem = ObstacleFlow(c)
    # solver.solve(problem)
    if len(sys.argv) == 4 and sys.argv[1] == "apps_config":
        hops_config()
    elif len(sys.argv) == 4:
        apps_routine(3*64, 64)
    else:
        c = nice_coeff(10)
        problem = ObstacleFlow(c,imax=3*64,jmax=64)
        solver.solve(problem)




if __name__ == '__main__':
    main()

