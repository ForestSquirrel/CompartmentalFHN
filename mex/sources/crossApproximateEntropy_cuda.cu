#include "mex.h"
#include <math.h>
#include "crossApproxEn.cu"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // Check number of inputs and outputs
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("crossApproxEntropy_cuda:invalidNumInputs", "Four inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("crossApproxEntropy_cuda:invalidNumOutputs", "One output required.");
    }

    // Parse inputs
    const mxArray* xSeries_mx = prhs[0];
    const mxArray* ySeries_mx = prhs[1];
    double* h_xSeries = mxGetPr(xSeries_mx);
    double* h_ySeries = mxGetPr(ySeries_mx);
    if (h_xSeries == nullptr || h_ySeries == nullptr) {
        mexErrMsgIdAndTxt("crossApproxEntropy_cuda:invalidInput", "Time series must be of type double.");
    }
    mwSize N_x = mxGetNumberOfElements(xSeries_mx);
    mwSize N_y = mxGetNumberOfElements(ySeries_mx);

    if (N_x != N_y) {
        mexErrMsgIdAndTxt("crossApproxEntropy_cuda:invalidInput", "Time series must be of the same length.");
    }
    mwSize N = N_x;

    int m = static_cast<int>(mxGetScalar(prhs[2]));
    double r_factor = mxGetScalar(prhs[3]);

    // Compute Cross-Approximate Entropy using the host function
    double crossApEn;
    crossApproxEnCUDA(h_xSeries, h_ySeries, N, m, r_factor, crossApEn);

    // Set the output
    plhs[0] = mxCreateDoubleScalar(crossApEn);
}
