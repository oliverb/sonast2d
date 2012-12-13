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

        /* border cells are always obstacles */
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
                FLAGS2(i,j) = (1-GEO2(i,j))*FC 
                            + (i==1)?:(1-GEO2(i-1,j))*B_W
                            + (i==imax)?:(1-GEO2(i+1,j))*B_E
                            + (j==1)?:(1-GEO2(i,j-1))*B_S
                            + (j==jmax)?:(1-GEO2(i,j+1))*B_N;
                flc += 1-GEO2(i,j);
            }
        }

        return_val = flc;
    """

    code2 = """
        #line 44 "flagfield.py"
        int imax = Ngeo[0]-2;
        int jmax = Ngeo[1]-2;

        /* border cells are always obstacles */
        for(int i=1; i<=imax; i++) {
            FLAGS2(i,0) = OC + (1-GEO2(i,1))*B_N;
            FLAGS2(i,jmax+1) = OC + (1-GEO2(i,jmax))*B_S;
        }

        for(int j=1; j<=jmax; j++) {
            FLAGS2(0,j) = OC + (1-GEO2(1,j))*B_E;
            FLAGS2(imax+1,j) = OC + (1-GEO2(imax,j))*B_W;
        }

        /* first pass */
        int flc = 0; // count fluid cells
        for(int i=1; i<=imax; i++) {
            for(int j=1; j<=jmax; j++) {
                if(GEO2(i,j)) {
                    FLAGS2(i,j) = OC;
                } else {
                    FLAGS2(i,j) = FC;
                    flc++;
                }
            }
        }

        /* second pass to update boundary information */
        int flc = 0; // count fluid cells
        for(int i=1; i<=imax; i++) {
            for(int j=1; j<=jmax; j++) {
                if(!(FLAG2(i,j)&FC)) {
                    FLAG2(i,j) += (FLAG2(i-1,j)&FC)*B_W
                                 +(FLAG2(i+1,j)&FC)*B_E
                                 +(FLAG2(i,j+1)&FC)*B_N
                                 +(FLAG2(i,j-1)&FC)*B_S;
                }
                /* TODO check flagfield for illegal obstacle cells */
            }
        }

        return_val = flc;
    """

    nfc = weave.inline(code2, ['flags', 'geo'], define_macros=cmacros.flag_dict)

    return nfc, flags
