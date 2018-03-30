#include "mex.h"
//#include "matrix.h"

#include <stdlib.h>
#include <stdio.h>
#include "surfpack_c_interface.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    /* variable declarations here */

    char *name;
    char *data_filename;
    int n_predictors;
    int n_responses;
    char *sps_filename;

    /* code here */

    name = mxArrayToString(prhs[0]);
    data_filename = mxArrayToString(prhs[1]);
    n_predictors = mxGetScalar(prhs[2]);
    n_responses = mxGetScalar(prhs[3]);
    sps_filename = mxArrayToString(prhs[4]);

    surfpack_load_data(name, data_filename, n_predictors, n_responses, 0);
    surfpack_build_model(name, "kriging");
    surfpack_save_model(name,sps_filename);

}


