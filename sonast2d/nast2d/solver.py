import numpy
import pylab

from nast2d.fd import time
from nast2d.fd import boundary
from nast2d.fd import pressure
from nast2d.fd import velocity
from nast2d.fd import flagfield

from nast2d.output import vtk

def solve(problem):
    xlength = problem.xlength
    ylength = problem.ylength
    Re = problem.Re
    geometry = problem.geo
    binfo = problem.boundary

    imax = problem.imax
    jmax = problem.jmax

    omega = problem.omega
    tau = problem.tau
    alpha = problem.alpha
    eps = problem.eps
    max_it = problem.max_it

    dx = xlength/imax
    dy = ylength/jmax


    ghost_dim = (imax+2, jmax+2)
    u = numpy.zeros(ghost_dim)
    v = numpy.zeros(ghost_dim)
    p = numpy.zeros(ghost_dim)
    rhs = numpy.zeros(ghost_dim)
    f = numpy.zeros(ghost_dim)
    g = numpy.zeros(ghost_dim)


    num_fc, flag = flagfield.generate_flagfield(geometry)
    print flag.shape

    t = 0.0
    iteration = 1
    acc_poisson_it = 0
    t_end = problem.t_end
    output_n = 0
    output_delta = 1.0

    Y = numpy.array([(j-1)*dy+dy/2.0 for j in range(1, jmax+1)])

    while t < t_end:
        dt = time.adaptive_stepwidth(dx, dy, u, v, Re, tau)
        boundary.set_outer_boundary(u, v, binfo)
        boundary.set_obstacle_slip(u, v, flag)
        velocity.velocity_guess(u, v, f, g, flag, dt, dx, dy, alpha, Re)
        velocity.compute_rhs(f, g, rhs, dt, dx, dy)
        it, res = pressure.solve_equation(p, rhs, flag, num_fc, dx, dy, eps, max_it, omega)
        velocity.velocity_correction(u, v, p, f, g, flag, dt, dx ,dy)

        if t >= output_n*output_delta:
            print "fuu"
            vtk.output_vector(u, v, flag, dx, dy, output_n, "test")
            #pylab.plot(Y, u[31,1:jmax+1], label=("t=%(time)f" % {"time":t}))
            output_n += 1

        print iteration, t, dt, it, res
        acc_poisson_it += it

        t += dt
        iteration += 1

#    for j in range(1, jmax+1):
#        print (j-1)*dy+dy/2.0, u[31, j]

    # Y = numpy.array([(j-1)*dy+dy/2.0 for j in range(1, jmax+1)])

    # ref_U = [1.0000, 0.8486, 0.7065, 0.5917, 0.5102, 
    #          0.4582, 0.4276, 0.4101, 0.3993, 0.3913, 
    #          0.3838, -0.0620, -0.3756, -0.3869, -0.3854, 
    #          -0.3690, -0.3381, -0.2960, -0.2472, -0.1951, 
    #          -0.1392, -0.0757, 0.0000] 
    # ref_U.reverse()
    # ref_U = numpy.array(ref_U)

    # ref_Y = [1.000, 0.990, 0.980, 0.970, 0.960, 0.950, 
    #          0.940, 0.930, 0.920, 0.910, 0.900, 0.500, 
    #          0.200, 0.180, 0.160, 0.140, 0.120, 0.100,
    #          0.080, 0.060, 0.040, 0.020, 0.000]
    # ref_Y.reverse()
    # ref_Y = numpy.array(ref_Y)

    # pylab.plot(Y, u[31,1:jmax+1], label="final")
    # pylab.plot(ref_Y, ref_U, label="reference")
    # pylab.axis('equal')
    # pylab.legend()
    # pylab.show()


