function freq_response = NFR(chirp_signal, system_response, fs, fft_opts)
% NFR computes the frequency response of a system using a chirp signal input.
% 
% Inputs:
%  - chirp_signal: The input chirp signal (1D array)
%  - system_response: The output signal from the system (1D array)
%  - fs: sampling rate
%  - fft_opts: Struct with optional fields for FFT analysis:
%       - 'window': Window function for STFT (default: hamming)
%       - 'nfft': Number of FFT points (default: length of window)
%       - 'noverlap': Overlap between segments (default: 50%)
%       - 'sides': 'onesided' or 'twosided' or 'centered' FFT (default: 'onesided')
%
% Output:
%  - freq_response: Magnitude and phase of system's frequency response
%    approximated from chirp input and system response

    % Default FFT options
    if nargin < 4
        fft_opts = struct();
    end
    window = getfield(fft_opts, 'window', hamming(1024));
    nfft = getfield(fft_opts, 'nfft', length(window));
    noverlap = getfield(fft_opts, 'noverlap', round(length(window) / 2));
    sides = getfield(fft_opts, 'sides', 'onesided');

    % Compute STFT of both chirp input and system response
    [S_chirp, F, T] = stft(chirp_signal, fs, 'Window', window, 'FFTLength', nfft, 'OverlapLength', noverlap, 'FrequencyRange', sides);
    [S_response, ~, ~] = stft(system_response, fs, 'Window', window, 'FFTLength', nfft, 'OverlapLength', noverlap, 'FrequencyRange', sides);

    % Get magnitude response (approximate frequency response)
    H = S_response ./ S_chirp; % Element-wise division to get the response

    % Get magnitude and phase
    magnitude = abs(H);
    phase = angle(H);

    % Output struct
    freq_response.magnitude = magnitude;
    freq_response.phase = phase;
    freq_response.frequency = F;
    freq_response.time = T;
end

function value = getfield(struct, field, default)
    if isfield(struct, field)
        value = struct.(field);
    else
        value = default;
    end
end