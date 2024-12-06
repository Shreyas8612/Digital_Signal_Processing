%% Time and Frequency domain Analysis
% Load the Signal
[sonar_signal, fs] = audioread('Sonar_Signal.wav');

% Inspect the outgoing Sonar Pulse and rest of the received signal
t = (0:length(sonar_signal)-1) / fs; % Time Vector

% Since the sampling frequency of the signal is 500 kHz
% For 1ms we sample at 10^-3 * 500kHz of frequency
sonar_signal_samples = round(1e-3 * fs); % Number of samples for the first 1ms

% Plot the outgoing pulse
figure;
plot(t(1:sonar_signal_samples), sonar_signal(1:sonar_signal_samples));
xlabel('Time (s)');
ylabel('Amplitude');
title('Outgoing Sonar Pulse in Time Domain');
grid on;

% Plot the entire received signal
figure;
plot(t, sonar_signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Received Sonar Signal in Time Domain');
grid on;

% Frequency Domain Analysis using pspectrum()
figure;
[pxx, f_ps] = pspectrum(sonar_signal(1:sonar_signal_samples), fs);
plot(f_ps / 1e3, 10*log10(pxx)); % Convert frequency to kHz and power to dB
xlabel('Frequency (kHz)');
ylabel('Power Spectrum (dB)');
title('Frequency Spectrum');
grid on;

% Frequency Domain Analysis using fft() fast fourier transforms
N = length(sonar_signal);
K = (0:1:N-1)';

% Apply Hamming window
H = 0.54 - (0.46 * (cos(2 * pi * (K / N))));
Windowed_Signal = sonar_signal .* H;

% Perform FFT on windowed and non-windowed signals
fft_signal_window = fft(Windowed_Signal);
fft_signal = fft(sonar_signal); % Compute FFT
f = (0:N-1) * (fs / N); % Frequency associated with each FFT bin.

% Normalization and logarithmic scaling
fft_magnitude = abs(fft_signal) / N; % Normalize FFT for non-windowed signal
fft_magnitude_window = abs(fft_signal_window) / N; % Normalize FFT for windowed signal

% Convert to dB for better representation
fft_magnitude_db = 20 * log10(fft_magnitude);
fft_magnitude_window_db = 20 * log10(fft_magnitude_window);

% Plot FFT magnitude for original and windowed signals
figure;
plot(f(1:N/2) / 1e3, fft_magnitude_db(1:N/2), 'b');
hold on;
plot(f(1:N/2) / 1e3, fft_magnitude_window_db(1:N/2), 'r');
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Frequency Domain Representation of Sonar Signal (Rectangular vs Hanning Window)');
legend('Rectangular Window', 'Hanning Window');
grid on;
hold off;

% Identify the carrier frequency
[~, peak_index] = max(fft_magnitude(1:N/2));
carrier_frequency = f(peak_index);
fc = carrier_frequency;
fprintf('The carrier frequency is %.1f kHz\n', carrier_frequency / 1e3);

%% Bandpass filter
fc_h = fc - 10000; % Cut-off frequency for High-pass filter
fc_l = fc + 10000; % Cut-off frequency for Low-pass filter

% Normalized Cut-Off Frequency
Omega_c_h = 2 * pi * fc_h / fs;
Omega_c_l = 2 * pi * fc_l / fs;

% Pre-warping the Cut-Off Frequency
w1 = tan(Omega_c_h / 2);
w2 = tan(Omega_c_l / 2);

% Define Analog Prototype of Butterworth Filter
syms s_ s z;
H_s = 1 / (s^2 + sqrt(2)*s + 1); % Analog prototype of 2nd order Butterworth filter

% Frequency Transformation for Band-pass Filter
H_s_hp = subs(H_s, s, (s_^2 + (w1 * w2)) / ((w2 - w1) * s_));

% Substitute Bilinear Transformation (s = (z-1)/(z+1))
s_bilinear = ((z - 1) / (z + 1));
H_z = subs(H_s_hp, s_, s_bilinear);

% Multiply numerator and denominator by (z+1)^2 to clear fractions
[num, den] = numden(H_z); % Get numerator and denominator of H(z)

% Convert to Coefficients for Standard IIR Form
[num_coeff, ~] = coeffs(num, z, 'All');
[den_coeff, ~] = coeffs(den, z, 'All');
b = double(num_coeff / den_coeff(1));
a = double(den_coeff / den_coeff(1));

% Ensure b and a are row vectors
b = reshape(b, 1, []);
a = reshape(a, 1, []);

% Butterworth Bandpass Filter Design using butter() function
cutoff_freqs = [(fc - 10e3), (fc + 10e3)]; % Bandpass cut-off frequencies in Hz

% Normalize the cut-off frequencies
Wn = cutoff_freqs / (fs / 2);

% Design a 2nd order Butterworth bandpass filter
[b_, a_] = butter(2, Wn, 'bandpass');

% Display the coefficients of butter() using Manual low-level code
disp('Manual computation of coefficients:');
disp('b (numerator):');
disp(b);
disp('a (denominator):');
disp(a);
disp('Manual computation of coefficients:');
disp('b (numerator):');
disp(b_);
disp('a (denominator):');
disp(a_);

% Implement Low-Level IIR Convolutional sum
Low_level_IIR_Output = zeros(size(Windowed_Signal)); % Initialize output

Z_1 = 0; % First delay element
Z_2 = 0; % Second delay element
Z_3 = 0; % Third delay element
Z_4 = 0; % Fourth delay element

Z_in_1 = 0; % First input delay element
Z_in_2 = 0; % Second input delay element
Z_in_3 = 0; % Third input delay element
Z_in_4 = 0; % Fourth input delay element

for n = 1:length(Windowed_Signal)
    % Calculate the current output sample using the canonical form realization
    Low_level_IIR_Output(n) = (b(1) * Windowed_Signal(n) + b(2) * Z_in_1 + b(3) * Z_in_2 + b(4) * Z_in_3 + b(5) * Z_in_4) ...
                                - (a(2) * Z_1 + a(3) * Z_2 + a(4) * Z_3 + a(5) * Z_4);
    % Update delay elements
    Z_in_4 = Z_in_3;
    Z_in_3 = Z_in_2;
    Z_in_2 = Z_in_1;
    Z_in_1 = Windowed_Signal(n);

    Z_4 = Z_3;
    Z_3 = Z_2;
    Z_2 = Z_1;
    Z_1 = Low_level_IIR_Output(n);
end

% Remove any NaN or Inf values that may have arisen
Low_level_IIR_Output(~isfinite(Low_level_IIR_Output)) = 0;

% Plot the Un-Filtered signal in time domain
figure;
t_wind = (0:length(Windowed_Signal)-1) / fs; % Time vector
plot(t_wind, Windowed_Signal, 'b', 'LineWidth', 0.5);
hold on;

% Plot the Filtered signal in time domain
plot(t_wind, Low_level_IIR_Output, 'r', 'LineWidth', 0.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Comparison of Filtered and Un-filtered signals in Time domain');
legend('Un-Filtered','Filtered', 'Location', 'Best');
grid on;
hold off;

% Plot the Filtered signal in frequency domain
figure;
[pxx_manual, f_ps_manual] = pspectrum(Low_level_IIR_Output, fs);
plot(f_ps_manual / 1e3, 10*log10(pxx_manual));
hold on;

% Plot the Un-Filtered signal in frequency domain
[pxx, f_ps] = pspectrum(Windowed_Signal, fs);
plot(f_ps / 1e3, 10*log10(pxx)); % Convert frequency to kHz and power to dB
xlabel('Frequency (kHz)');
ylabel('Power Spectrum (dB)');
title('Comparison of Filtered and Un-filtered signals in Frequency domain');
legend('Filtered Signal', 'Un-Filtered Signal', 'Location', 'Best');
grid on;

% FFT of filtered and windowed filtered signals
fft_filtered_signal_window = fft(Low_level_IIR_Output);

% Normalization for FFT magnitudes
fft_magnitude_windowed_filtered = abs(fft_filtered_signal_window) / N;

% Convert to dB for better representation
fft_magnitude_windowed_filtered_db = 20 * log10(fft_magnitude_windowed_filtered);

% Make sure N/2 is always an integer by rounding down
N_half = ceil(N / 2); % Use floor to ensure it's an integer

% Plot FFT magnitude for windowed original and windowed filtered signals
figure;
plot(f(1:N_half) / 1e3, fft_magnitude_window_db(1:N_half), 'b'); % Original windowed signal
hold on;
plot(f(1:N_half) / 1e3, fft_magnitude_windowed_filtered_db(1:N_half), 'r'); % Filtered windowed signal
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Hanning Window FFT: Original vs Filtered Signal');
legend('Windowed Original Signal', 'Windowed Filtered Signal');
grid on;
hold off;

%% Envelope detection
rectified_signal = abs(Low_level_IIR_Output);

% Plot Rectified Signal in Time Domain
figure;
t_rect = (0:length(rectified_signal)-1) / fs; % Time vector
plot(t_rect, rectified_signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Rectified Sonar Signal in Time Domain');
grid on;

% Plot Rectified Signal in Frequency Domain
figure;
[pxx_rectified, f_ps_rectified] = pspectrum(rectified_signal, fs);
plot(f_ps_rectified / 1e3, 10*log10(pxx_rectified)); % Convert frequency to kHz and power to dB
xlabel('Frequency (kHz)');
ylabel('Power Spectrum (dB)');
title('Rectified Sonar Signal in Frequency Domain');
grid on;

% Construction of an Ideal Impulse response
% Specifications
passband_edge = 1e3; % Passband edge frequency in Hz
stopband_edge = 2e3; % Stopband edge frequency in Hz
passband_ripple = 0.1; % Passband ripple in dB
stopband_attenuation = 50; % Stopband attenuation in dB

% Calculate Cutoff Frequency
f_cut = (passband_edge + stopband_edge) / 2; % Midpoint between passband and stopband edges
Fc = f_cut / fs; % Normalized cutoff frequency with respect to sampling frequency

% Calculate Number of Taps (N) using Hamming Window Parameters
% Transition Width = 3.3 / N
delta_f = stopband_edge - passband_edge; % Transition width in Hz
delta_f_normalized = delta_f / (fs / 2);
N = ceil(3.3 / delta_f_normalized); % Ceiling to ensure it's an integer

% Calculate M
M = (N - 1) / 2; % M is half of (N-1)

% Display calculated values
fprintf('Calculated Number of Taps (N): %d\n', N);
fprintf('Calculated Filter Order (M): %d\n', M);

% Generate Ideal Impulse Response of Lowpass Filter (Similar to Professor's Approach)
h_positive = zeros(1, M); % Preallocate for positive indices

for n = 1:M
    h_positive(n) = 2 * Fc * sin(n * 2 * pi * Fc) / (n * 2 * pi * Fc); % Ideal impulse response for LP filter (+ve n)
end

% Construct Full Impulse Response (Including n = 0 and Negative Half)
h_ideal = [fliplr(h_positive), 2 * Fc, h_positive]; % Add n = 0 coefficient and negative indices

% Plot Ideal Impulse Response
figure;
plot(h_ideal);
title('Ideal Impulse Response (Filter Coefficients)');
xlabel('FIR Filter Taps');
ylabel('Amplitude');
grid on;


% Generate Window Function (Hamming)
%K = (0:1:N-1)';
%window = 0.54 - (0.46 * (cos(2 * pi * (K / N)))); % Generate Hamming window
window = hamming(N)';
% Apply Window to Ideal Impulse Response
h_fir = h_ideal .* window; % Apply window to truncate the ideal impulse response

% Plot Windowed Filter Coefficients
figure;
plot(h_fir);
title('Windowed Filter Coefficients (Hamming Window Applied)');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% Low-Level FIR Convolution Implementation
y_fir_manual = zeros(size(rectified_signal)); % Initialize output
for q = 1:length(rectified_signal)
    % Convolution sum for FIR filter
    for k = 1:N
        if (q - k + 1) > 0
            y_fir_manual(q) = y_fir_manual(q) + h_fir(k) * rectified_signal(q - k + 1);
        end
    end
end

% Plot Manually Filtered Rectified Signal in Time Domain
figure;
plot(t_rect, y_fir_manual);
xlabel('Time (s)');
ylabel('Amplitude');
title('Manually Filtered Rectified Sonar Signal in Time Domain');
grid on;

% Plot Manually Filtered Rectified Signal in Frequency Domain
figure;
[pxx_manual_filtered, f_ps_manual_filtered] = pspectrum(y_fir_manual, fs);
plot(f_ps_manual_filtered / 1e3, 10*log10(pxx_manual_filtered));
xlabel('Frequency (kHz)');
ylabel('Power Spectrum (dB)');
title('Manually Filtered Rectified Sonar Signal in Frequency Domain');
grid on;

% Calculate Time Delay and Range of Target
% Convert the time threshold to samples to ignore the initial sonar pulse
ignore_samples = 1000; % Number of samples to ignore (equivalent to 1 ms at 500 kHz)

% Find the peak of the filtered rectified signal after the initial pulse
[~, max_index] = max(y_fir_manual(ignore_samples+1:end)); % Find peak after the initial 1 ms

% Calculate the Time Delay of the Echo Signal
% Time Delay = Time of Echo - Time of Initial Pulse (ignoring 1 ms)
time_delay = (max_index / fs) - (ignore_samples / fs); % Time delay in seconds

% Calculate Range to Target
speed_of_sound = 1500; % Speed of sound in water (m/s)
range = (time_delay * speed_of_sound) / 2; % Divide by 2 for two-way travel

% Display Results
fprintf('Time Delay: %.6f seconds\n', time_delay);
fprintf('Range to Target: %.2f meters\n', range);