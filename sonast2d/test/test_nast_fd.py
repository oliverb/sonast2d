import unittest

import numpy

from nast2d.fd import time
from nast2d.fd import pressure

class TestPressureSolver(unittest.TestCase):
    def setUp(self):
        self.imax = 2
        self.jmax = 2
        self.dx = 0.5
        self.dy = 0.5
        self.flag = 0x10 * numpy.ones( (self.imax+2, self.jmax+2), 
                                       dtype=numpy.int32)
        self.rhs = numpy.zeros( (self.imax+2, self.jmax+2))
        self.p = numpy.ones( (self.imax+2, self.jmax+2) )
        self.nfc = self.imax*self.jmax

    def test_pressure_norm(self):
        # Norm should be one, since it is normalized by fluid cells 
        pnorm = pressure._compute_pressure_norm(self.p, self.flag, self.nfc)
        self.assertAlmostEqual(pnorm, 1.0)

    def test_normalized_residual(self):
        # Residual should be zero, all constants are solutions for RHS zero
        pnorm = pressure._compute_pressure_norm(self.p, self.flag, self.nfc)

        residual = pressure._compute_normalized_residual(self.p, pnorm, 
                        self.rhs, self.flag, self.nfc, self.dx, self.dy);

        self.assertAlmostEqual(residual, 0.0)


    def test_zero_rhs(self):
        # Modify starting value to be non constant
        self.p[1,1] = 2.0
        self.p[2,2] = 2.0
        it, res = pressure.solve_equation(self.p, self.rhs, self.flag,
                    self.nfc, self.dx, self.dy, 1E-8, 800, 1.7)

        # Then test if we achieve a constant again
        non_ghost = self.p[1:self.imax+1, 1:self.jmax+1]
        diff = non_ghost - non_ghost[0,0]
        err = sum(sum(diff))
        self.assertAlmostEqual(err, 0.0)

    # TODO add test for rhs zero with simple geometry


class TestAdaptiveStepwidth(unittest.TestCase):
    def setUp(self):
        self.dx = 0.5
        self.dy = 0.5
        self.u = numpy.zeros( (4, 4) )
        self.v = numpy.zeros( (4, 4) )
        self.tau = 0.5
        self.Re = 1.0

    def test_minwidth(self):
        delta_t = time.adaptive_stepwidth(self.dx, self.dy, self.u, self.v,
                                          self.Re, self.tau)

        target = self.tau*self.Re/2.0/(1.0/self.dx/self.dx+1.0/self.dy/self.dy)

        self.assertAlmostEqual(delta_t, target)

if __name__ == '__main__':
    unittest.main()