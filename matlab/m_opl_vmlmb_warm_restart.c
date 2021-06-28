/*
 * m_opl_vmlmb_warm_restart.c -
 *
 * C-Matlab interface to OptimPack library function opl_vmlmb_warm_restart.
 */

#include <stdint.h>

#include "mex.h"
#include "optimpacklegacy.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1) {
        mexErrMsgTxt("expecting exactly 1 input argument");
    }
    if (nlhs != 1) {
        mexErrMsgTxt("expecting exactly 1 output");
    }
    opl_vmlmb_workspace_t* ws = (opl_vmlmb_workspace_t*)mxGetPr(prhs[0]);
    mwSize dims[] = {1, 1};
    plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    int32_t* task = (int32_t*)mxGetPr(plhs[0]);
    task[0] = opl_vmlmb_warm_restart(ws);
}
