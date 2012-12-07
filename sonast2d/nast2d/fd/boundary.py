from scipy import weave
from .. import cmacros

# For testing only no slip condition
def set_outer_boundary(u, v, flag):
    code = """
        #line 7 "boundary.py"
        int imax = Nu[0]-2;
        int jmax = Nu[1]-2;

        /* Left wall noslip */
        for(int j=1; j<=jmax; j++) {
            U2(0,j) = 0.0;
            V2(0,j) = -1.0*V2(1,j);
        }

        /* Left wall slip 
        for(int j=1; j<=jmax; j++) {
            U2(0,j) = 0.0;
            V2(0,j) = V2(1,j);
        } */

        /* Right wall noslip */
        for(int j=1; j<=jmax; j++) {
            U2(imax,j) = 0.0;
            V2(imax+1,j) = -1.0*V2(imax,j);
        }

        /* Right wall slip
        for(int j=1; j<=jmax; j++) {
            U2(imax,j) = 0.0;
            V2(imax+1,j) = V2(imax,j);
        } */

        /* Top wall noslip */
        for(int i=1; i<=imax; i++) {
            V2(i,jmax) = 0.0;
            U2(i,jmax+1) = -1.0*U2(i,jmax);
        }

        /* Top wall slip
        for(int i=1; i<=imax; i++) {
            V2(i,jmax) = 0.0;
            U2(i,jmax+1) = U2(i,jmax);
        } */

        /* Bottom wall noslip */
        for(int i=1; i<=imax; i++) {
            V2(i,0) = 0.0;
            U2(i,0) = -1.0*U2(i,1);
        }

        /* Bottom wall slip
        for(int i=1; i<=imax; i++) {
            V2(i,0) = 0.0;
            U2(i,0) = U2(i,1);
        } */

        /* Driven cavity */
        for(int i=1; i<=imax; i++) {
            U2(i,jmax+1) = 2.0-U2(i,jmax);
        }
        """ 
    weave.inline(code, ['u', 'v', 'flag'])

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
