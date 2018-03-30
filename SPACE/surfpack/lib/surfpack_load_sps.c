#include "mex.h"
//#include "matrix.h"

#include <stdlib.h>
#include <stdio.h>
#include "surfpack_c_interface.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    /* variable declarations here */

    char *name;
    char *sps_filename;

    /* code here */

    name = mxArrayToString(prhs[0]);
    sps_filename = mxArrayToString(prhs[1]);

    surfpack_load_model(name,sps_filename);

}


