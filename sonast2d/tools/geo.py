import sys

import numpy

from nast2d.output import vtk

def make_bitfield(chi, dim=2, imax=10, jmax=10, kmax=10, 
        xlength=1.0, ylength=1.0, zlength=1.0):
    """Make 0/1 bitmask from indicator function:
        mask[i,j]=chi(i*dx,j*dy)"""
    if dim == 2:
        mask = numpy.zeros((imax+2,jmax+2), dtype=numpy.int32)
    elif dim == 3:
        mask = numpy.zeros((imax+2,jmax+2,kmax+2), dtype=numpy.int32)
    else:
        print "dim only 2 or 3..."
        return

    dx, dy, dz = xlength/imax, ylength/jmax, zlength/kmax

    for i in range(1, imax+1):
        print "i=%d of %d" % (i, imax)
        for j in range(1, jmax+1):
            if dim == 2:
                mask[i,j] = chi(numpy.r_[(i-1.)*dx+.5*dx, (j-1.)*dy+.5*dy])
            elif dim == 3:
                for k in range(1, kmax+1):
                    mask[i,j,k] = chi(numpy.r_[(i-1.)*dx+.5*dx, (j-1.)*dy+.5*dy,
                                        (k-1.)*dz+.5*dz]) 

    return mask

def output_bitfield_3d(mask, dx, dy, dz, filename="default.vtk"):
    m1, m2, m3 = mask.shape
    imax = m1-2
    jmax = m2-2
    kmax = m3-2

    with open(filename, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Scalar Field\n")
        f.write("ASCII\n")
        f.write("DATASET RECTILINEAR_GRID\n")

        f.write("DIMENSIONS %d %d %d\n" % (imax, jmax, kmax, ))

        f.write("X_COORDINATES %d double\n" % (imax, ))
        for i in range(1, imax+1):
            f.write("%f " % (dx*(i-1.+.5)))
        f.write("\n")

        f.write("Y_COORDINATES %d double\n" % (jmax, ))
        for j in range(1, jmax+1):
            f.write("%f " % (dy*(j-1.+.5)))
        f.write("\n")

        f.write("Z_COORDINATES %d double\n" % (kmax, ))
        for k in range(1, kmax+1):
            f.write("%f " % (dz*(k-1.+.5)))
        f.write("\n")

        f.write("POINT_DATA %d\n" % (imax*jmax*kmax, ))
        f.write("Scalars Mask int\n")
        f.write("LOOKUP_TABLE default\n")

        for k in range(1, kmax+1):
            for j in range(1, jmax+1):
                for i in range(1, imax+1):
                    if mask[i,j,k]:
                        f.write("1\n")
                    else:
                        f.write("0\n")

# Ball in center of [0,1]^3 of radius .25 ...
def _test_chi(x):
    dim, = x.shape
    
    dist = 0.0
    for i in range(0, dim):
        dist += (.5-x[i])*(.5-x[i])

    if dist <= 0.25*0.25:
        return 1
    return 0

def main():
    test = make_bitfield(_test_chi, dim=3, imax=100, jmax=100, kmax=100)
    test = 0x10*test
    output_bitfield_3d(test, 1./10, 1./10, 1./10, filename="fooo.vtk")

if __name__ == '__main__':
    main()
