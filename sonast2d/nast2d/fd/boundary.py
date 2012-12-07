from scipy import weave
from .. import cmacros

# For testing only no slip condition
def set_outer_boundary(u, v, flag, boundary={}):
    # north = boundary.get('north', 'noslip')
    # east = boundary.get('east', 'noslip')
    # south = boundary.get('south', 'noslip')
    # west = boundary.get('west', 'noslip')

    m, n = u.shape
    imax = m - 2
    jmax = m - 2

    # Left wall noslip
    u[0,1:jmax+1] = 0.0
    v[0,1:jmax+1] = -1.0*v[1,1:jmax+1]

    # Right wall noslip
    u[imax,1:jmax+1] = 0.0
    v[imax+1,1:jmax+1] = -1.0*v[imax,1:jmax+1]

    # Top wall noslip
    v[1:imax+1,jmax] = 0.0
    u[1:imax+1,jmax+1] = -1.0*u[1:imax+1,jmax]

    # Bottom wall noslip
    v[1:imax+1,0] = 0.0
    u[1:imax+1,0] = -1.0*u[1:imax+1,1]

    # Driven cavity
    u[1:imax+1,jmax+1] = 2.0-u[1:imax+1,jmax]

def set_obstacle_slip(u, v, flag):
    code = """
        #line 44 "boundary.py"
        int imax = Nu[0]-2;
        int jmax = Nu[1]-2;

        for(int i=1; i<imax; i++) {
            for(int j=1; j<jmax; j++) {
                if( !(FLAG2(i,j)&FC) && FLAG2(i,j)!=OC ) {
                    switch(FLAG2(i,j)) {
                    case B_N:
                        V2(i,j) = 0.0;
                        U2(i,j) = -1.0*U2(i,j+1);
                        U2(i-1,j) = -1.0*U2(i-1,j+1);
                        break;
                    case B_E:
                        U2(i,j) = 0.0;
                        V2(i,j) = -1.0*V2(i+1,j);
                        V2(i,j-1) = -1.0*V2(i+1,j-1);
                        break;
                    case B_S:
                        V2(i,j-1) = 0.0;
                        U2(i,j) = -1.0*U2(i,j-1);
                        U2(i-1,j) = -1.0*U2(i-1,j-1);
                        break;
                    case B_W:
                        U2(i-1,j) = 0.0;
                        V2(i,j) = -1.0*V2(i-1,j);
                        V2(i,j-1) = 1.0*V2(i-1,j-1);
                        break;
                    case B_NE:
                        U2(i,j) = 0.0;
                        V2(i,j) = 0.0;
                        U2(i-1,j) = -1.0*U2(i-1,j+1);
                        V2(i,j-1) = -1.0*V2(i+1,j-1);
                        break;
                    case B_SE:
                        U2(i,j) = 0.0;
                        V2(i,j-1) = 0.0;
                        U2(i-1,j) = -1.0*U2(i-1,j-1);
                        V2(i,j) = -1.0*V2(i+1,j);
                        break;
                    case B_SW:
                        U2(i-1,j) = 0.0;
                        V2(i,j-1) = 0.0;
                        U2(i,j) = -1.0*U2(i,j-1);
                        V2(i,j) = -1.0*V2(i-1,j);
                        break;
                    case B_NW:
                        U2(i-1,j) = 0.0;
                        V2(i,j) = 0.0;
                        V2(i,j-1) = -1.0*V2(i-1,j-1);
                        U2(i,j) = -1.0*U2(i,j+1);
                        break;
                    }
                }
            }
        }
        """
    weave.inline(code, ['u', 'v', 'flag'], define_macros=cmacros.flag_dict)
