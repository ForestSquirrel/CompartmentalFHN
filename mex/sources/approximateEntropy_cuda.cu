#include "mex.h"
#include <math.h>
#include "approxEn.cu"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // Check number of inputs and outputs
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("approxEntropy_cuda:invalidNumInputs", "Three inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("approxEntropy_cuda:invalidNumOutputs", "One output required.");
    }

    // Parse inputs
    const mxArray* timeSeries_mx = prhs[0];
    double* h_timeSeries = mxGetPr(timeSeries_mx);
    if (h_timeSeries == nullptr) {
        mexErrMsgIdAndTxt("approxEntropy_cuda:invalidInput", "timeSeries must be of type double.");
    }
    mwSize N = mxGetNumberOfElements(timeSeries_mx);

    int m = static_cast<int>(mxGetScalar(prhs[1]));
    double r_factor = mxGetScalar(prhs[2]);

    // Compute Approximate Entropy using the host function
    double ApEn;
    approxEnCUDA(h_timeSeries, N, m, r_factor, ApEn);

    // Set the output
    plhs[0] = mxCreateDoubleScalar(ApEn);
}
