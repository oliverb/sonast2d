from .. import cmacros

def output_scalar(p, flags, dx, dy, n, name="default_p"):
    m1, m2 = p.shape
    imax = m1-2
    jmax = m2-2

    filename = "%s_%04d.vtk" % (name, n, )
    with open(filename, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Scalar Field\n")
        f.write("ASCII\n")
        f.write("DATASET RECTILINEAR_GRID\n")

        f.write("DIMENSIONS %d %d 1\n" % (imax, jmax, ))

        f.write("X_COORDINATES %d double\n" % (imax, ))
        for i in range(1, imax+1):
            f.write("%f " % (dx*(i-1)))
        f.write("\n")

        f.write("Y_COORDINATES %d double\n" % (jmax, ))
        for j in range(1, jmax+1):
            f.write("%f " % (dy*(j-1)))
        f.write("\n")

        f.write("Z_COORDINATES 1 double\n0.0\n")

        f.write("POINT_DATA %d\n" % (imax*jmax, ))
        f.write("Scalars Pressure double\n")
        f.write("LOOKUP_TABLE default\n")

        for j in range(1, jmax+1):
            for i in range(1, imax+1):
                if flags[i,j]&cmacros.FC:
                    f.write("%f\n" % (p[i,j]))
                else:
                    f.write("0.0\n")

def output_vector(u, v, flags, dx, dy, n, name="default_v"):
    #filename = "default_{0}".format(n)
    m1, m2 = u.shape
    imax = m1-2
    jmax = m2-2

    filename = "%s_%04d.vtk" % (name, n, )
    with open(filename, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Vector Field\n")
        f.write("ASCII\n")
        f.write("DATASET RECTILINEAR_GRID\n")

        f.write("DIMENSIONS %d %d 1\n" % (imax, jmax, ))

        f.write("X_COORDINATES %d double\n" % (imax, ))
        for i in range(1, imax+1):
            f.write("%f " % (dx*(i-1)))
        f.write("\n")

        f.write("Y_COORDINATES %d double\n" % (jmax, ))
        for j in range(1, jmax+1):
            f.write("%f " % (dy*(j-1)))
        f.write("\n")

        f.write("Z_COORDINATES 1 double\n0.0\n")

        f.write("POINT_DATA %d\n" % (imax*jmax, ))
        f.write("VECTORS Velocity double\n")

        for j in range(1, jmax+1):
            for i in range(1, imax+1):
                if flags[i,j]&cmacros.FC:
                    u_av = 0.5*(u[i-1, j]+u[i, j])
                    v_av = 0.5*(v[i, j-1]+v[i, j])
                    f.write("%f %f 0.0\n" % (u_av, v_av))
                else:
                    f.write("0.0 0.0 0.0\n")
