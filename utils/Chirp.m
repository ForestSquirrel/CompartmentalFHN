function [t, chirp_signal, f_t] = Chirp(A, tmax, fs, f0, f1, options)
    % Generates a chirp signal and an array of frequencies based on the input parameters
    %
    % Arguments:
    %   A      - Amplitude of the signal (positive scalar)
    %   tmax   - End time of the signal [s] (positive scalar)
    %   fs     - Sampling rate [Hz] (positive scalar)
    %   f0     - Starting frequency [Hz] (positive scalar)
    %   f1     - Final frequency [Hz] (positive scalar)
    %
    % Options:
    %   type   - Type of chirp ('log', 'linear', 'exp') (default: 'log')
    %   signal - Base function ('sin', 'cos', 'square') (default: 'sin')
    %
    % Output:
    %   t            - Time vector
    %   chirp_signal - Generated chirp signal array
    %   f_t          - Frequency at each time point array

    arguments
        A (1, 1) {mustBePositive}
        tmax (1, 1) {mustBePositive}
        fs (1, 1) {mustBePositive}
        f0 (1, 1) {mustBePositive}
        f1 (1, 1) {mustBePositive}
        options.type {mustBeMember(options.type, {'log','linear','exp'})} = 'log'
        options.signal {mustBeMember(options.signal, {'sin','cos','square'})} = 'sin'
    end

    % Generate time vector
    t = 0:1/fs:tmax;

    % Select chirp frequency variation
    if f0 == f1
        f_t = f0 * ones(size(t));
    else
        switch options.type
            case 'linear'
                f_t = f0 + (f1 - f0) * (t / tmax);
            case 'exp'
                f_t = f0 * (f1 / f0) .^ (t / tmax);
            case 'log'
                k = log(f1 / f0) / tmax;
                f_t = f0 * exp(k * t);
            otherwise
                error('Invalid chirp type');
        end
    end

    % Generate the base function
    phase = 2 * pi * cumtrapz(t, f_t);
    switch options.signal
        case 'sin'
            chirp_signal = A * sin(phase);
        case 'cos'
            chirp_signal = A * cos(phase);
        case 'square'
            chirp_signal = A * square(phase);
        otherwise
            error('Invalid signal type');
    end
end
