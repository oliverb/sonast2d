from .. import cmacros
from scipy import weave

def adaptive_stepwidth(dx, dy, u, v, Re, tau):
    # Compute greates absolute velocity in both directions
    code_max_velocity = """
        #line 7 "time.py"
        int imax = Nu[0]-2;
        int jmax = Nu[0]-2;

        double max_u = 0.0;
        double max_v = 0.0;

        for(int i=1; i<imax; i++) {
            for(int j=1; j<jmax; j++) {
                if(fabs(U2(i, j)) > max_u) {
                    max_u = fabs(U2(i, j));
                }
                if(fabs(V2(i, j)) > max_v) {
                    max_v = fabs(V2(i, j));
                }
            }
        }

        py::tuple results(2);
        results[0] = max_u;
        results[1] = max_v;
        return_val = results;
        """

    max_u, max_v = weave.inline(code_max_velocity, ['u', 'v'],
                                headers=['<cmath>'], libraries=['m'],
                                define_macros=cmacros.flag_dict);

    tmp = Re/2.0/(1.0/dx/dx+1.0/dy/dy);


    if max_u != 0.0 and tmp >= dx/max_u:
        tmp = dx/max_u

    if max_v != 0.0 and tmp >= dy/max_v:
        tmp = dy/max_v

    return tau*tmp



