#include "mex.h"
#include "MexInterface.h"

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[] )
{
    if(__mexFunction__(nlhs,plhs,nrhs,prhs) ==EXIT_FAILURE)
        mexErrMsgTxt("Fast Marching failure");
}