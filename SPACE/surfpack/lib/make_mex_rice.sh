#!/bin/sh

mex surfpack_load.c CFLAGS='-shared -rdynamic $CFLAGS' -I/home/juliendm/UTILS/Dakota/exe/include -L/home/juliendm/UTILS/Dakota/exe/lib -lsurfpack_c_interface -lsurfpack -lsurfpack_fortran -lncsuopt -lconmin -lboost_serialization -lpecos -lpecos_src -lteuchos -llapack -lblas -lstdc++
mex surfpack_eval.c CFLAGS='-shared -rdynamic $CFLAGS' -I/home/juliendm/UTILS/Dakota/exe/include -L/home/juliendm/UTILS/Dakota/exe/lib -lsurfpack_c_interface -lsurfpack -lsurfpack_fortran -lncsuopt -lconmin -lboost_serialization -lpecos -lpecos_src -lteuchos -llapack -lblas -lstdc++
mex surfpack_load_sps.c CFLAGS='-shared -rdynamic $CFLAGS' -I/home/juliendm/UTILS/Dakota/exe/include -L/home/juliendm/UTILS/Dakota/exe/lib -lsurfpack_c_interface -lsurfpack -lsurfpack_fortran -lncsuopt -lconmin -lboost_serialization -lpecos -lpecos_src -lteuchos -llapack -lblas -lstdc++
