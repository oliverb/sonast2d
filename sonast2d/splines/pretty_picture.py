from bspline import BasisSpline
import numpy
import pylab as plt

spline = BasisSpline(3, numpy.array([0.0, 4.0/3.0, 8.0/3.0, 4.0]))
ball = numpy.array([0.0, 0.8, 0.979795, 0.979795, 0.8, 0.0])

test_x = numpy.array(range(0,6))/5.0*4.0

# Basis values
mat = []
for i in range(0, 6):
    bval = spline.compute_basis_values(test_x[i])
    mat.append(numpy.copy(bval))
mat = numpy.array(mat)

circ = lambda x: numpy.sqrt(4.0-numpy.power(x-2.0, 2))
b = circ(test_x)

c = numpy.linalg.solve(mat, b)
ball = c

# Generate spline values
x_values = []
basis_values = []
ball_values = []
N = 100
dx = 4.0/N
for i in range(0, N+1):
    x = i*dx
    x_values.append(x)
    bval = spline.compute_basis_values(x)
    basis_values.append(numpy.copy(bval*ball))
    ball_values.append(spline.evaluate(ball,x))

basis_values = numpy.array(basis_values)
x_values = numpy.array(x_values)
ball_values = numpy.array(ball_values)

m, n = basis_values.shape
for i in range(0, n):
    plt.plot(x_values, basis_values[:,i])
plt.plot(x_values, ball_values)

plt.vlines(4.0/3.0, 0.0, circ(4.0/3.0), colors='k', linestyles='dashed')

plt.show()

