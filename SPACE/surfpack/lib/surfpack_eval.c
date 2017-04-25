#include "mex.h"

#include <stdlib.h>
#include <stdio.h>
#include "surfpack_c_interface.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    /* variable declarations */

    char *name;
    double *dvs;
    size_t ncols;
    double coeff;

    name = mxArrayToString(prhs[0]);
    dvs = mxGetPr(prhs[1]);
    ncols = mxGetN(prhs[1]);

    /* code */

    coeff = surfpack_eval_model(name, dvs, ncols);

    plhs[0] = mxCreateDoubleScalar(coeff);

}
