#ifndef MexInterface_h
#define MexInterface_h

#pragma GCC visibility push(default)

#include "mex.h"

EXTERN_C
int __mexFunction__(int nlhs, mxArray *plhs[],
                     int nrhs, const mxArray *prhs[] );
#pragma GCC visibility pop

#endif