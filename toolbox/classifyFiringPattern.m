function firingPattern = classifyFiringPattern(response, fs, MinPeakHeight, constants)
% Classifies the firing pattern of a neuron based on inter-spike intervals (ISI),
% first spike latency (fsl), post-spike silence (pss), and slow after hyperpolarizing wave amplitude (swa).
%
% Reference pseudocode:
% https://hippocampome.org/php/Help_Firing_Pattern_Identification_Pseudocode.php
%
% INPUT:
% response - response of a neuron
% fs - sampling rate of response
% MinPeakHeight - for peak finding
% constants (optional) - constants for algortihm
% OUTPUT:
% firingPattern - string describing the firing pattern

% Constants
if nargin < 4
    MIN_SWA = 5; % mV
    DELAY_FACTOR = 2;
    TSTUT_PRE_FACTOR = 2.5;
    TSTUT_POST_FACTOR = 1.5;
    PSTUT_FACTOR = 5;
    SLN_FACTOR = 2;
    SLOPE_THRESHOLD = 0.003;
    MIN_TSTUT_FREQ = 25; % Hz
else
    MIN_SWA = constants.MIN_SWA;
    DELAY_FACTOR = constants.DELAY_FACTOR;
    TSTUT_PRE_FACTOR = constants.TSTUT_PRE_FACTOR;
    TSTUT_POST_FACTOR = constants.TSTUT_POST_FACTOR;
    PSTUT_FACTOR = constants.PSTUT_FACTOR;
    SLN_FACTOR = constants.SLN_FACTOR;
    SLOPE_THRESHOLD = constants.SLOPE_THRESHOLD;
    MIN_TSTUT_FREQ = constants.MIN_TSTUT_FREQ;
end

% coefficient for discharge
discharge_coef = 7;
% factor for stutering
factor_stut = 1.5;

% Obtain response peaks
[~, locs] = findpeaks(response, fs, "MinPeakHeight", MinPeakHeight);

% ISI - array of inter-spike intervals
ISI = diff(locs);

if numel(ISI) < 5
    error("Not enough spikes to find pattern \n" + ...
        "Peak Finder found %d spikes \n" + ...
        "Minimum required numver of spikes to identify pattern is 6", numel(ISI)+1);
end

% ISI times
ISInTimes = locs(1:end-1)+0.5*ISI;

% fsl - first spike latency
firstNoNZeroResponseElement = find(response, 1, "first");
fsl = locs(1) - firstNoNZeroResponseElement;

% Check for slow wave inbetween spikes
if any(ISI(2:end-1) > factor_stut * ISI(1:end-2))
    % Get indices for the slow wave slice in between spikes
    largeISIIndex = find(ISI == max(ISI), 1);
    startIdx = locs(largeISIIndex);
    endIdx = locs(largeISIIndex + 1);
else
    % Slow wave after spikes
    startIdx = locs(end);
    endIdx = find(response, 1, "last")/fs;
end

% Call the separate function to compute swa and pss
[swa, pss] = computeSlowWave(response, startIdx, endIdx, fs, discharge_coef);

% Initialize pattern as empty
firingPattern = '';

% Step (i) Check for Delay (D.)
if hasDelay(ISI, fsl, DELAY_FACTOR)
    firingPattern = [firingPattern 'D.'];
end

% Step (ii) Check for Transient Stuttering (TSTUT.) or Slow-Wave Bursting (TSWB.)
if hasTSTUT(ISI, swa, pss, MIN_SWA, TSTUT_PRE_FACTOR, TSTUT_POST_FACTOR, MIN_TSTUT_FREQ)
    if swa > MIN_SWA
        firingPattern = [firingPattern 'TSWB.'];
    else
        firingPattern = [firingPattern 'TSTUT.'];
    end
end

% Step (iii) Run statistical tests for steady-state firing
[slope, sigDiff] = runSolverStatTests(ISI, ISInTimes);

if ~sigDiff(1) % No significant improvement from 1-parameter to 2-parameter fit
    firingPattern = [firingPattern 'NASP/STEADY_STATE.'];
    return;
elseif ~sigDiff(2) % No significant improvement from 2-parameter to 3-parameter fit
    if slope > SLOPE_THRESHOLD
        firingPattern = [firingPattern 'ASP.'];
    else
        firingPattern = [firingPattern 'NASP.'];
    end
    return;
elseif ~sigDiff(3) % No significant improvpement from 3-parameter to 4-parameter fit
    firingPattern = [firingPattern 'ASP.'];
    firingPattern = [firingPattern 'ASP.'];
    return;
end

% Step (iv) Check for Persistent Stuttering (PSTUT.) or Slow-Wave Bursting (PSWB.)
if hasPSTUT(ISI, PSTUT_FACTOR)
    if swa > MIN_SWA
        firingPattern = [firingPattern 'PSWB.'];
    else
        firingPattern = [firingPattern 'PSTUT.'];
    end
end

% Step (v) Check for Silence (SLN)
if hasSLN(pss, ISI, SLN_FACTOR)
    firingPattern = [firingPattern 'SLN.'];
end
end

%% Subfunctions:
function [swa, pss] = computeSlowWave(response, startIdx, endIdx, fs, discharge_coef)
    % Extract the slice of the response for slow wave analysis
    slowWaveSliceApprox = response((startIdx*fs):(endIdx*fs));
    
    % Compute post-spike silence (pss)
    slowWaveLastIndice = max( find(gradient(slowWaveSliceApprox) <= min(gradient(slowWaveSliceApprox))./discharge_coef, 1, "last"), ...
                              find(gradient(slowWaveSliceApprox) >= max(gradient(slowWaveSliceApprox))./discharge_coef, 1, "last"));
    pss = slowWaveLastIndice - startIdx;

    % Compute slow wave amplitude (swa)
    slowWaveSlice = response((startIdx*fs):slowWaveLastIndice);
    slowWaveMaximasIndices = islocalmax(slowWaveSlice);
    slowWaveMinimasIndices = islocalmin(slowWaveSlice);

    swa = ( ...
        abs(mean( slowWaveSlice(slowWaveMinimasIndices) )) + ...
        abs(mean( slowWaveSlice(slowWaveMaximasIndices) )) ...
        )/2;
end


function isDelay = hasDelay(ISI, fsl, DELAY_FACTOR)
% Checks if the firing pattern has a delay
isDelay = fsl > DELAY_FACTOR * mean(ISI(1:2));
end

function isTSTUT = hasTSTUT(ISI, swa, pss, MIN_SWA, TSTUT_PRE_FACTOR, TSTUT_POST_FACTOR, MIN_TSTUT_FREQ)
% Checks if the firing pattern has transient stuttering (TSTUT.)
isTSTUT = false;
for i = 2:4
    if ISI(i) > ISI(i-1) * TSTUT_PRE_FACTOR && ...
            ISI(i) > ISI(i+1) * TSTUT_POST_FACTOR && ...
            mean(ISI(i:end)) > mean(ISI(1:i-1)) * TSTUT_PRE_FACTOR && ...
            (1 / mean(ISI(1:i-1))) > MIN_TSTUT_FREQ

        isTSTUT = true;
        return;
    end
    if pss > ISI(end) * TSTUT_PRE_FACTOR && ...
            (1 / mean(ISI(1:end))) > MIN_TSTUT_FREQ && ...
            swa > MIN_SWA

        isTSTUT = true;
        return;
    end
end
end

function [slope, sigDiff] = runSolverStatTests(ISI, ISInTimes)
% Performs piecewise linear fits and returns the slope and significance tests

% Initialize sigDiff
sigDiff = [false false false];

% Define X and Y
X = ISInTimes(:); % Ensure column vector
Y = ISI(:);       % Ensure column vector

% Model 1: Y = c (constant model)
c1 = mean(Y);
Y_fit1 = c1 * ones(size(Y));
SSR1 = sum((Y - Y_fit1).^2);
DF1 = length(Y) - 1; % Degrees of freedom for model 1

% Model 2: Y = a*X + b (linear model)
p2 = polyfit(X, Y, 1);
a2 = p2(1);
b2 = p2(2);
Y_fit2 = polyval(p2, X);
SSR2 = sum((Y - Y_fit2).^2);
DF2 = length(Y) - 2; % Degrees of freedom for model 2

slope = a2; % Return the slope of the linear fit

% Model 3: Piecewise linear-constant model
% Y = m*X + n if X < (c - n)/m
% Y = c if X >= (c - n)/m

% Initial parameter estimates: [m, n, c]
params0 = [a2, b2, c1];

% Fit model 3 using nonlinear regression
[params3,~,~,~,~] = nlinfit(X, Y, @model3_func, params0);
Y_fit3 = model3_func(params3, X);
SSR3 = sum((Y - Y_fit3).^2);
DF3 = length(Y) - 3; % Degrees of freedom for model 3

% Model 4: Piecewise linear model with two segments
% Y = a1*X + b1 if X < (b2 - b1)/(a1 - a2)
% Y = a2*X + b2 if X >= (b2 - b1)/(a1 - a2)

% Initial parameter estimates: [a1, b1, a2, b2]
params0 = [a2, b2, a2, b2];

% Fit model 4 using nonlinear regression
[params4,~,~,~,~] = nlinfit(X, Y, @model4_func, params0);
Y_fit4 = model4_func(params4, X);
SSR4 = sum((Y - Y_fit4).^2);
DF4 = length(Y) - 4; % Degrees of freedom for model 4

% Perform F-tests between models
% Compare models 1 and 2
F12 = ((SSR1 - SSR2) / (DF1 - DF2)) / (SSR2 / DF2);
p12 = 1 - fcdf(F12, DF1 - DF2, DF2);

% Compare models 2 and 3
F23 = ((SSR2 - SSR3) / (DF2 - DF3)) / (SSR3 / DF3);
p23 = 1 - fcdf(F23, DF2 - DF3, DF3);

% Compare models 3 and 4
F34 = ((SSR3 - SSR4) / (DF3 - DF4)) / (SSR4 / DF4);
p34 = 1 - fcdf(F34, DF3 - DF4, DF4);

% Set significance thresholds
pValThresholds = [0.05; 0.025; 0.016]; % [1 & 2; 2 & 3; 3 & 4]

% Determine if improvements are significant
sigDiff(1) = p12 < pValThresholds(1);
sigDiff(2) = p23 < pValThresholds(2);
sigDiff(3) = p34 < pValThresholds(3);

% Nested functions for models 3 and 4
    function Y_pred = model3_func(params, X)
        m = params(1);
        n = params(2);
        c = params(3);
        X0 = (c - n) / m;
        Y_pred = zeros(size(X));
        idx = X < X0;
        Y_pred(idx) = m * X(idx) + n;
        Y_pred(~idx) = c;
    end

    function Y_pred = model4_func(params, X)
        a1 = params(1);
        b1 = params(2);
        a2 = params(3);
        b2 = params(4);
        X0 = (b2 - b1) / (a1 - a2);
        Y_pred = zeros(size(X));
        idx = X < X0;
        Y_pred(idx) = a1 * X(idx) + b1;
        Y_pred(~idx) = a2 * X(~idx) + b2;
    end
end


function isPSTUT = hasPSTUT(ISI, PSTUT_FACTOR)
% Checks for Persistent Stuttering (PSTUT.)
ISI_max = max(ISI);
factor_1 = ISI_max / ISI(find(ISI == ISI_max) - 1);
factor_2 = ISI_max / ISI(find(ISI == ISI_max) + 1);
isPSTUT = ...
    (factor_1 + factor_2) > PSTUT_FACTOR;
end

function isSLN = hasSLN(pss, ISI, SLN_FACTOR)
% Checks for Silence (SLN) after firing
isSLN = ...
    pss > SLN_FACTOR * mean(ISI(end-1:end)) && ...
    pss > SLN_FACTOR * max(ISI);
end
