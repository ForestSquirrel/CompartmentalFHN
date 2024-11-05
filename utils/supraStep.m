function [t, sig] = supraStep(Amplitude, Frequency, SamplingRate, NumPeriods, options)
% Generates a parametrized square wave signal
%
% Arguments:
%   Amplitude      - Amplitude of the wave (positive scalar)
%   Frequency      - Frequency of the wave [Hz] (positive scalar)
%   SamplingRate   - Sampling rate [Hz] (positive scalar)
%   NumPeriods     - Number of periods of the wave (positive integer)
%
% Options:
%   type                - Type of wave ('square', 'step', 'bipolar') [default: 'square']
%                         'square' : Normal square wave (-Amplitude to Amplitude)
%                         'step'   : Step wave (0 to Amplitude)
%                         'bipolar': Alternating wave (Amplitude, 0, -Amplitude, 0)
%   inverse             - Invert the signal (true/false) [default: false]
%   numZeroPeriodsStart - Number of zero periods at the start (non-negative integer) [default: 1]
%   numZeroPeriodsEnd   - Number of zero periods at the end (non-negative integer) [default: 1]
%
% Returns:
%   t   - Time vector
%   sig - Generated signal
%
% Example:
%   [t, sig] = GenerateWave(1, 1, 100, 2, 'type', 'bipolar', 'inverse', true);

    arguments
        Amplitude (1,1) {mustBePositive, mustBeFinite}
        Frequency (1,1) {mustBePositive, mustBeFinite}
        SamplingRate (1,1) {mustBePositive, mustBeFinite}
        NumPeriods (1,1) {mustBePositive, mustBeInteger, mustBeFinite}
        options.type {mustBeMember(options.type, {'square','step','bipolar'})} = 'square'
        options.inverse (1,1) logical = false
        options.numZeroPeriodsStart (1,1) {mustBeNonnegative, mustBeInteger, mustBeFinite} = 1
        options.numZeroPeriodsEnd (1,1) {mustBeNonnegative, mustBeInteger, mustBeFinite} = 1
    end

    % Compute period
    T = 1 / Frequency;

    % Compute total duration
    t_zero_start = options.numZeroPeriodsStart * T;
    t_signal = NumPeriods * T;
    t_zero_end = options.numZeroPeriodsEnd * T;
    t_total = t_zero_start + t_signal + t_zero_end;

    % Generate time vector
    dt = 1 / SamplingRate;
    t = 0:dt:t_total - dt; % Ensure the time vector does not exceed t_total

    % Initialize signal
    sig = zeros(size(t));

    % Indices for different sections
    N_zero_start = round(options.numZeroPeriodsStart * T * SamplingRate);
    N_signal = round(NumPeriods * T * SamplingRate);
    idx_signal = N_zero_start + 1 : N_zero_start + N_signal;

    % Time vector for the signal part
    t_sig = t(idx_signal) - t(idx_signal(1)); % Start from zero

    % Generate signal based on type
    switch options.type
        case 'square'
            sig(idx_signal) = Amplitude * square(2 * pi * Frequency * t_sig);
        case 'step'
            sig(idx_signal) = (Amplitude / 2) * (square(2 * pi * Frequency * t_sig, 50) + 1);
        case 'bipolar'
            t_frac = mod(t_sig, T) / T;
            sig_part = zeros(size(t_sig));
            idx1 = t_frac >= 0 & t_frac < 0.25;
            idx2 = t_frac >= 0.25 & t_frac < 0.5;
            idx3 = t_frac >= 0.5 & t_frac < 0.75;
            idx4 = t_frac >= 0.75 & t_frac < 1.0;
            sig_part(idx1) = Amplitude;
            sig_part(idx2) = 0;
            sig_part(idx3) = -Amplitude;
            sig_part(idx4) = 0;
            sig(idx_signal) = sig_part;
        otherwise
            error('Invalid type specified.');
    end

    % Apply inversion if specified
    if options.inverse
        sig = -sig;
    end
end
