import numpy
import pylab as py
from mpl_toolkits.mplot3d import Axes3D
from scipy import integrate

class BasisSpline:
    def __init__(self, order, knots):
        self.order = order
        num_knots, = knots.shape
        self.knots = numpy.zeros(num_knots+2*self.order)
        self.knots[0:self.order] = knots[0]
        self.knots[self.order+num_knots:] = knots[num_knots-1]
        self.knots[self.order:self.order+num_knots] = knots
        num_knots, = self.knots.shape
        self.num_knots = num_knots
        self.n = num_knots - order - 1
        self.data = numpy.zeros( (self.order+1, num_knots) ) 

    def compute_basis_values(self, x):
        for d in range(0, self.order+1):
            for i in range(0, self.n+self.order-d):
                if d==0:
                    if self.knots[i] <= x <= self.knots[i+1]:
                        self.data[d, i] = 1.0
                    else:
                        self.data[d, i] = 0.0
                else:
                    if abs(self.knots[i+d]-self.knots[i])>0.00001:
                        fct1 = (x-self.knots[i])/(self.knots[i+d]-self.knots[i])
                    else:
                        fct1 = 0.0
                    if abs(self.knots[i+1+d]-self.knots[i+1])>0.00001:
                        fct2 = (self.knots[i+1+d]-x)/(self.knots[i+1+d]-self.knots[i+1])
                    else:
                        fct2 = 0.0
                    self.data[d, i] = fct1*self.data[d-1, i]+fct2*self.data[d-1,i+1]
        return self.data[self.order,0:self.n]

    def compute_basis_matrix(self, xvals):
        B = [] 
        for x in xvals:
            # Need to make a copy, otherwise the b in B just gets
            # overwritten at the next iteration
            b = numpy.copy(self.compute_basis_values(x))
            B.append(b)

        return numpy.array(B)

    def compute_basis_integrals(self):
        I = numpy.zeros(self.n)
        for i in range(0, self.n):
            b_i = lambda x: self.compute_basis_values(x)[i]
            I[i], err = integrate.quad(b_i, self.knots[0], self.knots[self.num_knots-1])
        return I


    def evaluate(self, c, x):
        self.compute_basis_values(x)
        return sum(self.data[self.order,0:self.n]*c)

    def print_in_range(self):
        c = numpy.array([0.0, 0.8, 0.979795, 0.979795, 0.8, 0.0])

        num_knots, = self.knots.shape
        x = self.knots[0]

        delta_x = (self.knots[num_knots-1]-self.knots[0])/500.0

        while x <= self.knots[num_knots-1]:
            line = str(x) + " "
            foo = self.compute_basis_values(x)
            for i, b in enumerate(foo):
                line = line + str(b*c[i]) + " "
            line = line + str(sum(foo*c))
            print line
            x += delta_x
        # foo = self.compute_basis_values(x)
        # line = str(x) + " "
        # for i, b in enumerate(foo):
        #     line = line + str(b*c[i]) + " "
        # line = line + str(sum(foo*c))
        # print line

class BasisSplineSurface:
    def __init__(self, order, xknots, yknots):
        self.order = order
        self.xknots = xknots
        self.yknots = yknots

        self.xspline = BasisSpline(self.order, self.xknots)
        self.yspline = BasisSpline(self.order, self.yknots)

    def evaluate(self, a, p):
        if p.shape[0] != 2:
            print "Error, wrong dimension ..."
            return 0.0

        x = p[0]
        y = p[1]

        # Could redo this as tensor product ...
        # So much wasted cpu time :-(
        ev = 0.0
        xvals = self.xspline.compute_basis_values(x)
        for i, Nx in enumerate(xvals):
            # Sum of Y basis functions at y
            SY = self.yspline.evaluate(a[i,:], y)
            ev += Nx*SY

        return ev

if __name__ == '__main__':
    xk = numpy.array([0.0,  1.0])
    yk = numpy.array([0.0,  1.0])
    a = numpy.array([[0.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], 
                     [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]])

    sf = BasisSplineSurface(3, xk, yk)

    X = numpy.linspace(0.0, 1.0, 25)
    X, Y = numpy.meshgrid(X, X)

    Z = []
    for xs, ys in zip(X, Y):
        zs = []
        for x, y in zip(xs, ys):
            zs.append(sf.evaluate(a, numpy.array([x,y])))
        Z.append(numpy.array(zs))
    Z = numpy.array(Z)

    fig = py.figure()
    ax = Axes3D(fig)

    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='hot')

    py.show()






