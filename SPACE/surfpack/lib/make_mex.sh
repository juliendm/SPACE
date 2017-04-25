#!/bin/sh

mex surfpack_load.c CFLAGS='-shared -rdynamic $CFLAGS' -I/ADL/juliendm/UTILS/Dakota/exe/include -L/usr/lib64/atlas -L/ADL/juliendm/UTILS/boost_1_49_0/exe/lib -L/ADL/juliendm/UTILS/Dakota/exe/lib -lsurfpack_c_interface -lsurfpack -lsurfpack_fortran -lncsuopt -lconmin -lboost_serialization -lpecos -lteuchos -llapack -lblas -lstdc++
mex surfpack_eval.c CFLAGS='-shared -rdynamic $CFLAGS' -I/ADL/juliendm/UTILS/Dakota/exe/include -L/usr/lib64/atlas -L/ADL/juliendm/UTILS/boost_1_49_0/exe/lib -L/ADL/juliendm/UTILS/Dakota/exe/lib -lsurfpack_c_interface -lsurfpack -lsurfpack_fortran -lncsuopt -lconmin -lboost_serialization -lpecos -lteuchos -llapack -lblas -lstdc++
