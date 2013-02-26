from math import pi
import sys

import numpy
import scipy.optimize

from splines import bspline
from tools import geo

def obstacle_chi(c):
    N, = c.shape

    # For N degrees of freedom we need N-2 knots
    curve = bspline.BasisSpline(3, numpy.linspace(0.0, 4.0, N-2))

    def chi(x):
        if x[0] <= 2.0:
            return 0
        if x[0] >= 6.0:
            return 0

        sval = curve.evaluate(c, x[0]-2.0)

        print x[0], x[1], 4.0-sval, 4.0+sval

        if abs(x[1] - 4.0) <= sval:
            return 1
        else:
            return 0

    return chi

# Routines to calculate nice circulat starting values ...
def half_circle(x):
    if x >= 2.0:
        return 0.0
    x = x-1.0
    return numpy.sqrt(1.0-x*x)

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

def main():
    c = nice_coeff(6)
    print c
    chi = obstacle_chi(c)
    test = geo.make_bitfield(chi, dim=3, imax=3*64, jmax=64, kmax=64,
                            xlength=24.0, ylength=8.0, zlength=8.0)
    test = 0x10*test
    geo.output_bitfield_3d(test, 24./3./64., 8./64., 8./64., filename="fooo.vtk")

if __name__ == '__main__':
    main()
