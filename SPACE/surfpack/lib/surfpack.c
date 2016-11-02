/* Example of loading a Surfpack model from a surface file (.sps or
   .bsps) and evaluating it at a number of points. */

#include <stdlib.h>
#include <stdio.h>
#include "surfpack_c_interface.h"

void load_data(char *name, char *data_filename, int n_predictors, int n_responses, int n_cols_to_skip) {
	surfpack_load_data(name, data_filename, n_predictors, n_responses, n_cols_to_skip);
}

void save_data(char *name, char *data_filename) {
	surfpack_save_data(name, data_filename);
}

void load_model(char *name, char *model_filename) {
	surfpack_load_model(name, model_filename);
}

void save_model(char *name, char *model_filename) {
	surfpack_save_model(name, model_filename);
}

void add(char *name, double *point, int size, double f) {
	surfpack_add_data(name, point, size, f);
}

void build(char *name, char *type) {
	surfpack_build_model(name, type);
}

double eval(char *name, double *point, int size) {
    return surfpack_eval_model(name, point, size);
}

double variance(char *name, double *point, int size) {
    return surfpack_variance_model(name, point, size);
}
