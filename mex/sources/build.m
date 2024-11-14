%% BUILD mex 

%% Sample Entropy
%mexcuda -output sampleEntropy_cuda sampEnt_cuda.cu -I"D:\CompartmentalFHN\CUDA"

%% Approximate Entropy
%mexcuda -output approximateEntropy_cuda approximateEntropy_cuda.cu -I"D:\CompartmentalFHN\CUDA"

%% Cross-Approximate Entropy
%mexcuda -output crossApproximateEntropy_cuda crossApproximateEntropy_cuda.cu -I"D:\CompartmentalFHN\CUDA"