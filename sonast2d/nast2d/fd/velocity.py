from scipy import weave
from .. import cmacros

def velocity_guess(u, v, f, g, flag, dt, dx, dy, alpha, Re):
    code = """
        #line 6 "velocity_guess"
        int imax = Nu[0]-2;
        int jmax = Nu[1]-2;

        double dd_u_x = 0.0;
        double dd_u_y = 0.0;
        double d_u2_x = 0.0;
        double d_uv_y = 0.0;

        /* Inner points of F, (12) */
        for(int i=1; i<=imax-1; i++) {
            for(int j=1; j<=jmax; j++) {
                if( FLAG2(i, j)&FC && FLAG2(i+1, j)&FC ) {
                    /* D^2 u / D x^2  --- (4) */
                    dd_u_x = (U2(i+1, j) - 2.0* U2(i, j) + U2(i-1, j))/dx/dx;

                    /* D^2 u / D y^2  --- (4) */
                    dd_u_y = (U2(i, j+1) - 2.0* U2(i, j) + U2(i, j-1))/dy/dy;

                    /* D u^2 / D x  --- (4) */
                    d_u2_x = ((U2(i, j)+U2(i+1, j))*(U2(i,j)+U2(i+1,j))/4.0
                             -(U2(i-1,j)+U2(i, j))*(U2(i-1,j)+U2(i,j))/4.0)/dx;
                    d_u2_x += (fabs(U2(i,j)+U2(i+1,j))*(U2(i,j)-U2(i+1,j))/4.0
                              -fabs(U2(i-1,j)+U2(i,j))*(U2(i-1,j)-U2(i,j))/4.0)
                              *alpha/dx;

                    /* D uv / D y  --- (4) */
                    d_uv_y = ((V2(i,j)+V2(i+1,j))*(U2(i,j)+U2(i,j+1))/4.0
                            -(V2(i,j-1)+V2(i+1,j-1))*(U2(i,j-1)+U2(i,j))/4.0)/dy;
                    d_uv_y += (fabs(V2(i,j)+V2(i+1,j))*(U2(i,j)-U2(i,j+1))/4.0
                            -fabs(V2(i,j-1)+V2(i+1,j-1))*(U2(i,j-1)-U2(i,j))/4.0)
                            *alpha/dy;

                    F2(i,j) = U2(i,j) + dt*( (dd_u_x + dd_u_y)/Re - d_u2_x - d_uv_y );
                }
                // else {
                // F2(i,j) = U2(i,j);
                //}
            }
        }

        /* Inner points of G, (13) */
        for(int i=1; i<=imax; i++) {
            for(int j=1; j<=jmax-1; j++) {
                if( FLAG2(i,j)&FC && FLAG2(i,j+1)&FC ) {
                    /* D^2 v / D x^2 (5) */
                    dd_u_x = (V2(i+1, j) - 2.0*V2(i, j) + V2(i-1, j))/dx/dx;

                    /* D^2 v / D y^2 (5) */
                    dd_u_y = (V2(i, j+1) - 2.0*V2(i, j) + V2(i, j-1))/dy/dy;

                    /* D v^2 / D y (5) */
                    d_u2_x = ((V2(i,j)+V2(i,j+1))*(V2(i,j)+V2(i,j+1))/4.0
                            -(V2(i,j-1)+V2(i,j))*(V2(i,j-1)+V2(i,j))/4.0)/dy;  
                    d_u2_x += (fabs(V2(i,j)+V2(i,j+1))*(V2(i,j)-V2(i,j+1))/4.0
                              -fabs(V2(i,j-1)+V2(i,j))*(V2(i,j-1)-V2(i,j))/4.0)
                              *alpha/dy;

                    /* D uv / D x (5) */
                    d_uv_y = ((U2(i,j)+U2(i,j+1))*(V2(i,j)+V2(i+1,j))/4.0
                            -(U2(i-1,j)+U2(i-1,j+1))*(V2(i-1,j)+V2(i,j))/4.0)/dx;
                    d_uv_y += (fabs(U2(i,j)+U2(i,j+1))*(V2(i,j)-V2(i+1,j))/4.0
                              -fabs(U2(i-1,j)+U2(i-1,j+1))*(V2(i-1,j)-V2(i,j))/4.0)
                              *alpha/dx;

                    G2(i,j) = V2(i,j) + dt*( (dd_u_x + dd_u_y)/Re - d_u2_x - d_uv_y );
                }
                //  else {
                //      G2(i,j) = V2(i,j);
                // }
            }
        }

        /* Outer boundary (really neccessary?) */
        for(int j=1; j<=jmax; j++) {
            F2(0,j) = U2(0,j);
            F2(imax,j) = U2(imax,j);
        }
        for(int i=1; i<=imax; i++) {
            G2(i,0) = V2(i,0);
            G2(i,jmax) = V2(i,jmax);
        }

        /* Obstacle boundary */
        for(int i=1; i<=imax; i++) {
            for(int j=1; j<=jmax; j++) {
                if( !(FLAG2(i,j)&FC) && FLAG2(i,j)!=OC ) {
                    switch(FLAG2(i,j)) {
                    case B_N:
                        G2(i,j) = V2(i,j);
                        break;
                    case B_W:
                        F2(i-1,j) = U2(i-1,j);
                        break;
                    case B_E:
                        F2(i,j) = U2(i,j);
                        break;
                    case B_S:
                        G2(i,j-1) = V2(i,j-1);
                        break;
                    case B_NW:
                        G2(i,j) = V2(i,j);
                        F2(i-1,j) = U2(i-1,j);
                        break;
                    case B_SW:
                        F2(i-1,j) = U2(i-1,j);
                        G2(i,j-1) = V2(i,j-1);
                        break;
                    case B_SE:
                        G2(i,j-1) = V2(i,j-1);
                        F2(i,j) = U2(i,j);
                        break;
                    case B_NE:
                        G2(i,j) = V2(i,j);
                        F2(i,j) = U2(i,j);
                        break;    
                    }
                }
            }
        }
        """

    variables = "u v f g flag dt dx dy alpha Re".split(" ")
    weave.inline(code, variables, define_macros=cmacros.flag_dict,
                 headers=['<cmath>'], libraries=['m'])

def velocity_correction(u, v, p, f, g, flag, dt, dx, dy):
    code = """
        #line 133 "velocity.py"
        int imax = Nu[0]-2;
        int jmax = Nu[1]-2;

        for(int i=1; i<=imax-1; i++) {
            for(int j=1; j<=jmax; j++) {
                if( FLAG2(i, j)&FC && FLAG2(i+1, j)&FC ) {
                    U2(i, j) = F2(i, j) - dt/dx*(P2(i+1, j) - P2(i, j));
                }
            }
        }

        for(int i=1; i<=imax; i++) {
            for(int j=1; j<=jmax-1; j++) {
                if( FLAG2(i, j)&FC && FLAG2(i, j+1)&FC ) {
                    V2(i, j) = G2(i, j) - dt/dy*(P2(i, j+1) - P2(i, j));
                }
            }
        }
        """

    variables = "u v p f g flag dt dx dy".split(" ")
    weave.inline(code, variables, define_macros=cmacros.flag_dict)

def compute_rhs(f, g, rhs, dt, dx, dy):
    code = """
        #line 159 "velocity.py"
        int imax = Nf[0]-2;
        int jmax = Nf[1]-2;

        for(int i=1; i<=imax; i++) {
            for(int j=1; j<=jmax; j++) {
                RHS2(i, j) = (F2(i, j) - F2(i-1, j))/dx
                            +(G2(i, j) - G2(i, j-1))/dy;
                RHS2(i, j) = RHS2(i, j)/dt;
            }
        }
        """

    variables = "f g rhs dt dx dy".split(" ")
    weave.inline(code, variables, define_macros=cmacros.flag_dict)