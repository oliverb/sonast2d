import numpy
from scipy import weave

from .. import cmacros

def generate_flagfield(geo):
    m, n = geo.shape

    flags = numpy.zeros((m,n), dtype=numpy.int32)

    code = """
        #line 12 "flagfield.py"
        int imax = Ngeo[0]-2;
        int jmax = Ngeo[1]-2;

        /* border cells are always obstacles, but this is mainly
         * cosmetics ...
         */
        for(int i=1; i<=imax; i++) {
            FLAGS2(i,0) = OC + (1-GEO2(i,1))*B_N;
            FLAGS2(i,jmax+1) = OC + (1-GEO2(i,jmax))*B_S;
        }

        for(int j=1; j<=jmax; j++) {
            FLAGS2(0,j) = OC + (1-GEO2(1,j))*B_E;
            FLAGS2(imax+1,j) = OC + (1-GEO2(imax,j))*B_W;
        }

        /* inner cells */
        int flc = 0; // count fluid cells
        for(int i=1; i<=imax; i++) {
            for(int j=1; j<=jmax; j++) {
                if(GEO2(i,j) == 0) {
                    FLAGS2(i, j) = FC;
                    flc += 1;
                } else {
                    /* check for neighboring fluid cells to aid computation
                     * of boundary values
                     */
                    FLAGS2(i,j) = OC;
                    if(GEO2(i-1,j) == 0) {
                        FLAGS2(i,j) += B_W;
                    }
                    if(GEO2(i+1,j) == 0) {
                        FLAGS2(i,j) += B_E;
                    }
                    if(GEO2(i,j-1) == 0) {
                        FLAGS2(i,j) += B_S;
                    }
                    if(GEO2(i,j+1) == 0) {
                        FLAGS2(i,j) += B_N;
                    }
                }
            }
        }

        int ff_illegal=1;
        while(ff_illegal) {
            ff_illegal=0;
            /* "correct" flagfield */
            for(int i=1; i<=imax; i++) {
                for(int j=1; j<=jmax; j++) {
                    if(!(FLAGS2(i,j)&FC)) {
                        if(FLAGS2(i-1,j)&FC && FLAGS2(i+1,j)&FC ||
                           FLAGS2(i,j-1)&FC && FLAGS2(i,j+1)&FC) {
                            FLAGS2(i,j) += FC;
                            ff_illegal = 1;
                            flc += 1;
                        }
                    }
                }
            }           
        }


        return_val = flc;
    """

    nfc = weave.inline(code, ['flags', 'geo'], define_macros=cmacros.flag_dict)

    return nfc, flags
