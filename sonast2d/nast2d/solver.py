import numpy

from nast2d.fd import time
from nast2d.fd import boundary
from nast2d.fd import pressure
from nast2d.fd import velocity
from nast2d.fd import flagfield

from nast2d.fd import measurements

from nast2d.output import vtk

def solve(problem):
    # Ugly copy from problem class
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

        # At this point we compute the dissipation at time t, NOT t+dt,
        # since measurements.dissipation assumes that the obstacle cells
        # are set for enforcement of boundary conditions.
        dissipation = measurements.dissipation(u, v, flag, dx, dy)

        velocity.velocity_guess(u, v, f, g, flag, dt, dx, dy, alpha, Re)
        velocity.compute_rhs(f, g, rhs, dt, dx, dy)
        it, res = pressure.solve_equation(p, rhs, flag, num_fc, dx, dy, eps, max_it, omega)
        velocity.velocity_correction(u, v, p, f, g, flag, dt, dx ,dy)

        if t >= output_n*output_delta:
            vtk.output_vector(u, v, flag, dx, dy, output_n, "test")
            output_n += 1

        print iteration, t, dt, it, res, dissipation

        acc_poisson_it += it
        t += dt
        iteration += 1


