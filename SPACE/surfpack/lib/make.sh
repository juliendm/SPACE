#!/bin/sh

gcc surfpack.c -o surfpack.so -shared -fPIC -rdynamic -I ~/UTILS/Dakota/exe/include -L /usr/lib64/atlas -L ~/UTILS/boost_1_49_0/exe/lib -L ~/UTILS/Dakota/exe/lib -lsurfpack_c_interface -lsurfpack -lsurfpack_fortran -lncsuopt -lconmin -lboost_serialization -lpecos -lteuchos -llapack -lblas -lstdc++ -Wl,-rpath,~/UTILS/boost_1_49_0/exe/lib:~/UTILS/Dakota/exe/lib 

