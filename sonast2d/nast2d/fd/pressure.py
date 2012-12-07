from scipy import weave
from .. import cmacros

def solve_equation(p, rhs, flag, num_fc, dx, dy, eps, max_it, omega):
    iteration = 0
    residual = 1E10

    while iteration<max_it and residual>eps:
        if iteration%25 == 0:
            pressure_norm = _compute_pressure_norm(p, flag, num_fc)
            
        _set_boundary_conditions(p, flag, dx, dy)
        _relax(p, rhs, flag, dx, dy, omega)

        if iteration%25 == 0:
            residual = _compute_normalized_residual(p, pressure_norm, rhs, 
                                                    flag, num_fc, dx, dy)

        iteration += 1

    return (iteration, residual)

def _set_boundary_conditions(p, flag, dx, dy):
    code = """
        #line 19 "pressure.py"
        int imax = Np[0]-2;
        int jmax = Np[1]-2;

        /* Outer boundary */
        for(int i=1; i<=imax; i++) {
            P2(i, 0) = P2(i, 1);
            P2(i, jmax+1) = P2(i, jmax);
        }

        for(int j=1; j<=jmax; j++) {
            P2(0, j) = P2(1, j);
            P2(imax+1, j) = P2(imax, j);
        }

        /* Boundary conditions on obstacle cells */
        double tmp = 1.0/(dx*dx+dy*dy);

        for(int i=1; i<=imax; i++) {
            for(int j=1; j<=jmax; j++) {
                if( !(FLAG2(i, j)&FC) && !(FLAG2(i, j)!=OC)) {
                    switch(FLAG2(i, j)) {
                    case B_N:
                        P2(i, j) = P2(i, j+1);
                        break;
                    case B_E:
                        P2(i, j) = P2(i+1, j);
                        break;
                    case B_S:
                        P2(i, j) = P2(i, j-1);
                        break;
                    case B_W:
                        P2(i, j) = P2(i-1, j);
                        break;
                    case B_NE:
                        P2(i, j) = tmp*(dx*dx*P2(i, j+1)+dy*dy*P2(i+1, j));
                        break;
                    case B_SE:
                        P2(i, j) = tmp*(dx*dx*P2(i, j-1)+dy*dy*P2(i+1, j));
                        break;
                    case B_SW:
                        P2(i, j) = tmp*(dx*dx*P2(i, j-1)+dy*dy*P2(i-1, j));
                        break;
                    case B_NW:
                        P2(i, j) = tmp*(dx*dx*P2(i, j+1)+dy*dy*P2(i-1, j));
                        break;
                    }
                }
            }
        }
        """

    weave.inline(code, ['p', 'flag', 'dx', 'dy'],
                     define_macros=cmacros.flag_dict);

def _relax(p, rhs, flag, dx, dy, omega=1.7):
    code = """
        #line 76 "pressure.py"
        int imax = Np[0]-2;
        int jmax = Np[1]-2;

        double dx_inv_sq = 1.0/dx/dx;
        double dy_inv_sq = 1.0/dy/dy;
        double sor_factor = omega/2.0/(dx_inv_sq + dy_inv_sq);

        for(int i=1; i<=imax; i++) {
            for(int j=1; j<=jmax; j++) {
                if(FLAG2(i,j)&FC)
                P2(i, j) = (1.0 - omega)*P2(i, j) 
                         + sor_factor*( (P2(i+1, j)+P2(i-1, j))*dx_inv_sq
                                       +(P2(i, j+1)+P2(i, j-1))*dy_inv_sq
                                       -RHS2(i, j) );
            }
        }
        """

    variables = "p rhs flag dx dy omega".split(" ")

    weave.inline(code, variables, define_macros=cmacros.flag_dict)

def _compute_normalized_residual(p, pnorm, rhs, flag, num_fc, dx, dy):
    code = """
        #line 101 "pressure.py"
        int imax = Np[0]-2;
        int jmax = Np[1]-2;

        double dx_inv_sq = 1.0/dx/dx;
        double dy_inv_sq = 1.0/dy/dy;

        double tmp = 0.0;
        double residual = 0.0;

        for(int i=1; i<=imax; i++) {
            for(int j=1; j<=jmax; j++) {
                if(FLAG2(i, j)&FC)
                tmp = (P2(i+1, j)-2*P2(i, j)+P2(i-1, j))*dx_inv_sq
                     +(P2(i, j+1)-2*P2(i, j)+P2(i, j-1))*dy_inv_sq
                     -RHS2(i, j); 

                residual += tmp*tmp;
            }
        }

        residual = sqrt(residual/num_fc)/pnorm;

        return_val = residual;
        """

    variables = "p pnorm rhs flag num_fc dx dy".split(" ")

    return weave.inline(code, variables, define_macros=cmacros.flag_dict,
                        headers=['<cmath>'], libraries=['m'])


def _compute_pressure_norm(p, flag, num_fc):
    code = """
        #line 133 "pressure.py"
        int imax = Np[0]-2;
        int jmax = Np[1]-2;

        double acc = 0.0;

        for(int i=1; i<=imax; i++) {
            for(int j=1; j<=jmax; j++) {
                if(FLAG2(i, j)&FC)
                    acc += P2(i, j)*P2(i, j);
            }
        }

        double pnorm = sqrt(acc/num_fc);

        if(pnorm < 0.0001) {
            return_val = 1.0;
        } else {
            return_val = pnorm;
        }
        """

    return weave.inline(code, ['p', 'flag', 'num_fc'],
                        define_macros=cmacros.flag_dict,
                        headers=['<cmath>'], libraries=['m'])


