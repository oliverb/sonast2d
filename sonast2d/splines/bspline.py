import numpy
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

if __name__ == '__main__':
    test = BasisSpline(3, numpy.array([0.0, 4.0/3.0, 8.0/3.0, 4.0]))
    c = numpy.array([0.0, 0.8, 0.979795, 0.979795, 0.8, 0.0])
    bI = test.compute_basis_integrals()
    print sum(c*bI)
    func = lambda x: test.evaluate(c, x)
    cI, err = integrate.quad(func, 0.0, 4.0) 
    print cI




