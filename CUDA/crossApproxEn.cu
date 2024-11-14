#pragma once

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include "mex.h"

// CUDA error checking macro
#define cudaCheckError(ans) { gpuAssert((ans), __FILE__, __LINE__); }

/* SOURCE */
inline void gpuAssert(cudaError_t code, const char* file, int line, bool abort = true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "CUDA Error: %s %s %d\n",
            cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

__device__ inline bool isMatchXY(const double* xSeries, const double* ySeries, int i, int j, int m, double r)
{
    for (int k = 0; k < m; ++k)
    {
        if (fabs(xSeries[i + k] - ySeries[j + k]) > r)
            return false;
    }
    return true;
}

__global__ void computeCiXY(const double* xSeries, const double* ySeries, int N, int m, double r, double* d_Ci)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N - m + 1) return;

    int count = 0;
    int total = N - m + 1;

    for (int j = 0; j < total; ++j)
    {
        if (isMatchXY(xSeries, ySeries, i, j, m, r))
            count++;
    }

    d_Ci[i] = (double)count / total;
}

__host__ void crossApproxEnCUDA(const double* h_xSeries, const double* h_ySeries, int N, int m, double r_factor, double& crossApEn) {
    // Compute the mean and std_dev of both series concatenated
    double mean = 0.0;
    for (int i = 0; i < N; ++i) {
        mean += h_xSeries[i];
        mean += h_ySeries[i];
    }
    mean /= (2 * N);

    // Compute the standard deviation
    double std_dev = 0.0;
    for (int i = 0; i < N; ++i) {
        std_dev += (h_xSeries[i] - mean) * (h_xSeries[i] - mean);
        std_dev += (h_ySeries[i] - mean) * (h_ySeries[i] - mean);
    }
    std_dev = sqrt(std_dev / (2 * N));

    // Calculate tolerance `r` using the standard deviation
    double r = r_factor * std_dev;

    // Allocate device memory
    double* d_xSeries;
    double* d_ySeries;
    double* d_Ci_m;
    double* d_Ci_m1;
    cudaCheckError(cudaMalloc((void**)&d_xSeries, N * sizeof(double)));
    cudaCheckError(cudaMalloc((void**)&d_ySeries, N * sizeof(double)));
    cudaCheckError(cudaMalloc((void**)&d_Ci_m, (N - m + 1) * sizeof(double)));
    cudaCheckError(cudaMalloc((void**)&d_Ci_m1, (N - m) * sizeof(double)));

    // Copy data to device
    cudaCheckError(cudaMemcpy(d_xSeries, h_xSeries, N * sizeof(double), cudaMemcpyHostToDevice));
    cudaCheckError(cudaMemcpy(d_ySeries, h_ySeries, N * sizeof(double), cudaMemcpyHostToDevice));

    // Kernel launch parameters
    int threadsPerBlock = 256;
    int blocksPerGrid_m = (N - m + 1 + threadsPerBlock - 1) / threadsPerBlock;
    int blocksPerGrid_m1 = (N - m + threadsPerBlock - 1) / threadsPerBlock;

    // Launch kernels
    computeCiXY<<<blocksPerGrid_m, threadsPerBlock>>>(d_xSeries, d_ySeries, N, m, r, d_Ci_m);
    cudaCheckError(cudaGetLastError());

    computeCiXY<<<blocksPerGrid_m1, threadsPerBlock>>>(d_xSeries, d_ySeries, N, m + 1, r, d_Ci_m1);
    cudaCheckError(cudaGetLastError());

    // Copy results back to host
    double* h_Ci_m = new double[N - m + 1];
    double* h_Ci_m1 = new double[N - m];
    cudaCheckError(cudaMemcpy(h_Ci_m, d_Ci_m, (N - m + 1) * sizeof(double), cudaMemcpyDeviceToHost));
    cudaCheckError(cudaMemcpy(h_Ci_m1, d_Ci_m1, (N - m) * sizeof(double), cudaMemcpyDeviceToHost));

    // Compute Φ^xy(m, r) and Φ^xy(m + 1, r)
    double phi_m = 0.0;
    double phi_m1 = 0.0;

    for (int i = 0; i < N - m + 1; ++i)
    {
        if (h_Ci_m[i] > 0)
            phi_m += log(h_Ci_m[i]);
        else
            phi_m += log(1e-10); // Avoid log(0)
    }
    phi_m /= (N - m + 1);

    for (int i = 0; i < N - m; ++i)
    {
        if (h_Ci_m1[i] > 0)
            phi_m1 += log(h_Ci_m1[i]);
        else
            phi_m1 += log(1e-10); // Avoid log(0)
    }
    phi_m1 /= (N - m);

    // Compute cross-approximate entropy
    crossApEn = phi_m - phi_m1;

    // Free memory
    delete[] h_Ci_m;
    delete[] h_Ci_m1;
    cudaFree(d_xSeries);
    cudaFree(d_ySeries);
    cudaFree(d_Ci_m);
    cudaFree(d_Ci_m1);
}
