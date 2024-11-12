%% BUILD mex 
% mexcuda -output sampleEntropy_cuda sampEnt_cuda.cu -I"D:\CompartmentalFHN\CUDA"

t = 0:1e-4:5;

sine_wave = 1 * sin(2 * pi * t * 4);

tol = std(sine_wave) * 0.2;

tic
% https://www.mathworks.com/matlabcentral/fileexchange/35784-sample-entropy
m_sampEn = SampEn(2, tol, sine_wave)
toc

tic
cu_sampEn = sampleEntropy_cuda(sine_wave, 2, 0.2)
toc