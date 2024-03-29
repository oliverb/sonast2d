from scipy import weave
from .. import cmacros

# For testing only no slip condition
def set_outer_boundary(u, v, boundary={}):
    north = boundary.get('north', 'noslip')
    east = boundary.get('east', 'noslip')
    south = boundary.get('south', 'noslip')
    west = boundary.get('west', 'noslip')

    m, n = u.shape
    imax = m - 2
    jmax = n - 2

    if west == 'noslip':
        u[0,0:jmax+2] = 0.0
        v[0,0:jmax+2] = -1.0*v[1,0:jmax+2]
    elif west == 'slip':
        u[0,0:jmax+2] = 0.0
        v[0,0:jmax+2] = v[1,0:jmax+2]

    if east == 'noslip':
        u[imax,0:jmax+2] = 0.0
        v[imax+1,0:jmax+2] = -1.0*v[imax,0:jmax+2]
    elif east == 'slip':
        u[imax,0:jmax+2] = 0.0
        v[imax+1,0:jmax+2] = v[imax,0:jmax+2]

    if north == 'noslip':
        v[0:imax+2,jmax] = 0.0
        u[0:imax+2,jmax+1] = -1.0*u[0:imax+2,jmax]
    elif north == 'slip':
        v[0:imax+2,jmax] = 0.0
        u[0:imax+2,jmax+1] = u[0:imax+2,jmax]

    if south == 'noslip':
        v[0:imax+2,0] = 0.0
        u[0:imax+2,0] = -1.0*u[0:imax+2,1]
    elif south == 'slip':
        v[0:imax+2,0] = 0.0
        u[0:imax+2,0] = u[0:imax+2,1] 

    # Scenario specific boundary conditions
    if 'scenario' in boundary:
        if boundary['scenario'] == 'cavity':
            u[0:imax+2,jmax+1] = 2.0-u[0:imax+2,jmax]
        if boundary['scenario'] == 'obstacle':
            u[0,0:jmax+2] = 1.0
            v[0,0:jmax+2] = -1.0*v[1,0:jmax+2]

            u[imax,0:jmax+2] = 1.0
            v[imax,0:jmax+2] = -1.0*v[imax-1,0:jmax+2]

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
