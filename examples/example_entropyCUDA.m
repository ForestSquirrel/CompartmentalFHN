t = 0:1e-4:1;

sine_wave = 1 * sin(2 * pi * t * 4);

tol = std(sine_wave) * 0.2;

% Sample
tic
% https://www.mathworks.com/matlabcentral/fileexchange/35784-sample-entropy
m_sampEn = SampEn(2, tol, sine_wave)
toc

tic
cu_sampEn = sampleEntropy_cuda(sine_wave, 2, 0.2)
toc

% Approximate
tic
m_appEn = approximateEntropy(sine_wave)
%% 
toc

tic
cu_appEn = approximateEntropy_cuda(sine_wave, 2, 0.2)
toc

% Cross-Approximate
tic
cu_cross_appEn = crossApproximateEntropy_cuda(sine_wave, sine_wave, 2, 0.2)
toc

approximateEntropy_cuda()