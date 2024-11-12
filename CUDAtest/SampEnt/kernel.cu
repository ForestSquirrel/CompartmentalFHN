#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <math.h>

// CUDA error checking macro
#define cudaCheckError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char* file, int line, bool abort = true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "CUDA Error: %s %s %d\n",
            cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

// Device function to compute the distance between two vectors
__device__ inline bool isMatch(const double* timeSeries, int i, int j, int m, double r)
{
    for (int k = 0; k < m; ++k)
    {
        if (fabs(timeSeries[i + k] - timeSeries[j + k]) > r)
            return false;
    }
    return true;
}

// Kernel to compute counts for patterns of length m
__global__ void countMatches(const double* timeSeries, int N, int m, double r, unsigned long long* d_B)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N - m + 1) return;

    unsigned long long count = 0;
    for (int j = 0; j < N - m + 1; ++j)
    {
        if (i == j) continue; // Exclude self-matches
        if (isMatch(timeSeries, i, j, m, r))
            count++;
    }

    atomicAdd(d_B, count);
}

// Kernel to compute counts for patterns of length m+1
__global__ void countMatchesExtended(const double* timeSeries, int N, int m, double r, unsigned long long* d_A)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N - m) return;

    unsigned long long count = 0;
    for (int j = 0; j < N - m; ++j)
    {
        if (i == j) continue; // Exclude self-matches
        if (isMatch(timeSeries, i, j, m + 1, r))
            count++;
    }

    atomicAdd(d_A, count);
}

int main()
{
    // Input parameters
    const int N = 10000;        // Length of the time series
    const int m = 2;            // Embedding dimension
    const double r_factor = 0.2; // Tolerance as a fraction of standard deviation

    // Allocate and initialize host memory
    double* h_timeSeries = new double[N];
    // Initialize your time series data here
    for (int i = 0; i < N; ++i)
    {
        h_timeSeries[i] = sin(0.01 * i); // Example data
    }

    // Compute standard deviation for tolerance calculation
    double mean = 0.0;
    for (int i = 0; i < N; ++i) mean += h_timeSeries[i];
    mean /= N;

    double std_dev = 0.0;
    for (int i = 0; i < N; ++i) std_dev += (h_timeSeries[i] - mean) * (h_timeSeries[i] - mean);
    std_dev = sqrt(std_dev / N);

    double r = r_factor * std_dev;

    // Allocate device memory
    double* d_timeSeries;
    unsigned long long* d_A;
    unsigned long long* d_B;
    cudaCheckError(cudaMalloc((void**)&d_timeSeries, N * sizeof(double)));
    cudaCheckError(cudaMalloc((void**)&d_A, sizeof(unsigned long long)));
    cudaCheckError(cudaMalloc((void**)&d_B, sizeof(unsigned long long)));

    // Copy data to device
    cudaCheckError(cudaMemcpy(d_timeSeries, h_timeSeries, N * sizeof(double), cudaMemcpyHostToDevice));
    cudaCheckError(cudaMemset(d_A, 0, sizeof(unsigned long long)));
    cudaCheckError(cudaMemset(d_B, 0, sizeof(unsigned long long)));

    // Launch kernels
    int threadsPerBlock = 256;
    int blocksPerGrid_m = (N - m + 1 + threadsPerBlock - 1) / threadsPerBlock;
    int blocksPerGrid_m1 = (N - m + threadsPerBlock - 1) / threadsPerBlock;

    countMatches << <blocksPerGrid_m, threadsPerBlock >> > (d_timeSeries, N, m, r, d_B);
    cudaCheckError(cudaGetLastError());

    countMatchesExtended << <blocksPerGrid_m1, threadsPerBlock >> > (d_timeSeries, N, m, r, d_A);
    cudaCheckError(cudaGetLastError());

    // Copy results back to host
    unsigned long long h_A = 0;
    unsigned long long h_B = 0;
    cudaCheckError(cudaMemcpy(&h_A, d_A, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
    cudaCheckError(cudaMemcpy(&h_B, d_B, sizeof(unsigned long long), cudaMemcpyDeviceToHost));

    // Compute Sample Entropy
    if (h_B == 0 || h_A == 0)
    {
        printf("Sample Entropy is undefined (no matches found).\n");
    }
    else
    {
        double sampEn = -log((double)h_A / h_B);
        printf("Sample Entropy: %f\n", sampEn);
    }

    // Free device and host memory
    cudaFree(d_timeSeries);
    cudaFree(d_A);
    cudaFree(d_B);
    delete[] h_timeSeries;

    return 0;
}
