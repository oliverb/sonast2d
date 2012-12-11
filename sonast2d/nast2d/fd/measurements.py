from scipy import weave

from .. import cmacros

def dissipation(u, v, flag, dx, dy):
    # Computes dissipation at cell centers
    # C code computes at grid line intersections
    code = """
        #line 9 "measurements.py"
        int imax = Nu[0]-2;
        int jmax = Nu[0]-2;

        double dudx = 0.0;
        double dudy = 0.0;
        double dvdx = 0.0;
        double dvdy = 0.0;
        double u_t = 0.0;
        double u_b = 0.0;
        double v_l = 0.0;
        double v_r = 0.0;

        double dxi = 1.0/dx;
        double dyi = 1.0/dy;

        double dissipation_l2_sq = 0.0; 

        for(int i=1; i<=imax; i++) {
            for(int j=1; j<=jmax; j++) {
                if(FLAG2(i,j)&FC) {
                    dudx = dxi*(U2(i,j)-U2(i-1,j));
                    dvdy = dyi*(V2(i,j)-V2(i,j-1)); 

                    /* THE PLAN
                     * Bilinear interpolation of u to points at the middle
                     * top and bottom of the cell (i,j).
                     * This is where v normaly resides. If to the top or 
                     * bottom is an obstacle cell, use boundary value instead.
                     * 
                     * THE BUT
                     * Is this even a good idea? We then need all kinds of
                     * extra information. Bilinear Interpolation should
                     * interpolate to the boundary value anyway, thats
                     * the whole point in how dirichlet conditions are
                     * enforced in this code.
                     *
                     * MISSING POINT
                     * What the heck are SLIP conditions anyway? Maybe
                     * they need extra treatment?
                     */
                    u_t = 0.25*(U2(i-1,j+1)+U2(i-1,j)+U2(i,j+1)+U2(i,j));
                    u_b = 0.25*(U2(i-1,j-1)+U2(i-1,j)+U2(i,j-1)+U2(i,j));

                    dudy = dyi*(u_t-u_b);

                    /* Bilinear interpolation of v to left and right */
                    v_r = 0.25*(V2(i+1,j-1)+V2(i+1,j)+V2(i,j-1)+V2(i,j));
                    v_l = 0.25*(V2(i-1,j-1)+V2(i-1,j)+V2(i,j-1)+V2(i,j));

                    dvdx = dxi*(v_r-v_l);

                    dissipation_l2_sq += dx*dy*(dudx*dudx + dudy*dudy 
                                               +dvdx*dvdx + dvdy*dvdy);
                }
            }
        }

        return_val = dissipation_l2_sq;
        """

    return weave.inline(code, ['u', 'v', 'flag', 'dx', 'dy'], 
                 define_macros=cmacros.flag_dict)

