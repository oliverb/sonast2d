import numpy

from nast2d.fd import time
from nast2d.fd import boundary
from nast2d.fd import pressure
from nast2d.fd import velocity

from nast2d.output import vtk


def solve():
    xlength = 1.0
    ylength = 1.0
    imax = 64
    jmax = 64

    Re = 100.0
    omega = 1.7
    tau = 0.5
    alpha = 0.9
    eps = 1E-7
    max_it = 5000

    dx = xlength/imax
    dy = ylength/jmax

    ghost_dim = (imax+2, jmax+2)
    u = numpy.zeros(ghost_dim)
    v = numpy.zeros(ghost_dim)
    p = numpy.zeros(ghost_dim)
    rhs = numpy.zeros(ghost_dim)
    f = numpy.zeros(ghost_dim)
    g = numpy.zeros(ghost_dim)
    flag = 0x10 * numpy.ones(ghost_dim, dtype=numpy.int32)

    num_fc = imax*jmax


    t = 0.0
    iteration = 1
    t_end = 5.0
    output_n = 0
    output_delta = 0.04

    while t < t_end:
        dt = time.adaptive_stepwidth(dx, dy, u, v, Re, tau)
        boundary.set_outer_boundary(u, v, flag)
        boundary.set_obstacle_slip(u, v, flag)
        velocity.velocity_guess(u, v, f, g, flag, dt, dx, dy, alpha, Re)
        velocity.compute_rhs(f, g, rhs, dt, dx, dy)
        it, res = pressure.solve_equation(p, rhs, flag, num_fc, dx, dy, eps, max_it, omega)
        velocity.velocity_correction(u, v, p, f, g, flag, dt, dx ,dy)

        if t >= output_n*output_delta:
            vtk.output_vector(u, v, flag, dx, dy, output_n, "test")
            output_n += 1

        print iteration, t, dt, it, res

        t += dt
        iteration += 1

    for j in range(0, jmax+1):
        print j*dy, u[31, j]

