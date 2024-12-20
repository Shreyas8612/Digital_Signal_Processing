%% Time and frequency domain analysis

% Load the Signal
[sonar_signal, fs] = audioread('Sonar_Signal.wav');

% Inspect the outgoing Sonar Pulse and rest of the received signal
t_signal = (0:length(sonar_signal)-1) / fs; % Time Vector

% Since the sampling frequency of the signal is 500 kHz
% For 1ms we sample at 10^-3 * 500kHz of frequency
sonar_signal_samples = round(1e-3 * fs); % Number of samples for the first 1ms

% Plot the outgoing pulse
figure;
plot(t_signal(1:sonar_signal_samples), sonar_signal(1:sonar_signal_samples));
xlabel('Time (s)');
ylabel('Amplitude');
title('Outgoing Sonar Pulse in Time Domain');
grid on;

% Plot the entire received signal
figure;
plot(t_signal, sonar_signal);
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
N_window = 524288;
K = (0:1:N-1)';

% Apply Hamming or hanning window
%H = (0.5 * (1-cos(2*pi*K/(N))));
H = (0.54 - (0.46 * (cos(2 * pi * (K /N))))); % It retians more of the main frequency
Windowed_Signal = sonar_signal.* H;

% Perform FFT on windowed and non-windowed signals
fft_signal_window = fft(Windowed_Signal);
fft_signal = fft(sonar_signal); % Compute FFT
f = (0:N_window-1) * (fs / N_window); % Frequency associated with each FFT bin.

% Normalization and logarithmic scaling
fft_magnitude = abs(fft_signal) / N_window; % Normalize FFT for non-windowed signal
fft_magnitude_window = abs(fft_signal_window) / N_window; % Normalize FFT for windowed signal

% Convert to dB for better representation
fft_magnitude_db = 20 * log10(fft_magnitude);
fft_magnitude_window_db = 20 * log10(fft_magnitude_window);

% Plot FFT magnitude for original and windowed signals
figure;
plot(f(1:N_window/2) / 1e3, fft_magnitude_db(1:N_window/2), 'b');
hold on;
plot(f(1:N_window/2) / 1e3, fft_magnitude_window_db(1:N_window/2), 'r');
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('FFT Frequency Resolution = 1 (Rectangular vs Hanning Window)');
legend('Rectangular Window', 'Hamming Window');
grid on;
hold off;

% Find the frequency corresponding to the maximum power in the spectrum
[~, peak_index] = max(pxx); % Index of maximum power
carrier_frequency = f_ps(peak_index); % Carrier frequency in Hz
fc = carrier_frequency;

% Print the carrier frequency in kHz
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
H_s_hp = subs(H_s, s, (s_^2 + (w1 * w2)) / ((w2 - w1) * s_)); % Substitute s -> omega_ac / s to transform to high-pass

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

% Display the coefficients of butter() and Manual low-level code
disp('Manual computation of coefficients:');
disp('b (numerator):');
disp(b);
disp('a (denominator):');
disp(a);

disp('Coefficients from butter() function:');
disp('b (numerator):');
disp(b_);
disp('a (denominator):');
disp(a_);

% Uncomment the line below to use the built-in `filter()` function
%Low_level_IIR_Output = filter(b, a, sonar_signal);

Signal = sonar_signal;

% Implement Low-Level IIR Convolutional sum
Low_level_IIR_Output = zeros(size(Signal)); % Initialize output

Z_1 = 0; % First delay element
Z_2 = 0; % Second delay element
Z_3 = 0; % Third delay element
Z_4 = 0; % Fourth delay element

Z_in_1 = 0; % First input delay element
Z_in_2 = 0; % Second input delay element
Z_in_3 = 0; % Third input delay element
Z_in_4 = 0; % Fourth input delay element

for n = 1:length(Signal)
    % Calculate the current output sample using the canonical form realization
    Low_level_IIR_Output(n) = (b(1) * Signal(n) + b(2) * Z_in_1 + b(3) * Z_in_2 + b(4) * Z_in_3 + b(5) * Z_in_4) ...
                  - (a(2) * Z_1 + a(3) * Z_2 + a(4) * Z_3 + a(5) * Z_4);

    % Update delay elements
    Z_in_4 = Z_in_3;
    Z_in_3 = Z_in_2;
    Z_in_2 = Z_in_1;
    Z_in_1 = Signal(n);

    Z_4 = Z_3;
    Z_3 = Z_2;
    Z_2 = Z_1;
    Z_1 = Low_level_IIR_Output(n);
end

% Remove any NaN or Inf values that may have arisen
Low_level_IIR_Output(~isfinite(Low_level_IIR_Output)) = 0;

% Plot the Un-Filtered signal in time domain
figure;
t_wind = (0:length(Signal)-1) / fs; % Time vector
plot(t_wind, Signal, 'b', 'LineWidth', 0.5);
hold on;

% Plot the Filtered signal in time domain
plot(t_wind, Low_level_IIR_Output, 'r', 'LineWidth', 0.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Comparison of Filtered and Un-filtered signals in Time domain');
legend('Un-Filtered','Filtered', 'Location', 'Best');
grid on;
hold off;

% Visualize Frequency Response
[H_BP, W_BP] = freqz(Low_level_IIR_Output, 1, 4096, fs);
figure;
plot(W_BP / 1e3, 20*log10(abs(H_BP)));
xlabel('Frequency (kHz)');
ylabel('Magnitude');
title('Frequency Response of IIR Bandpass Filter');
grid on;

% Plot the Filtered signal in frequency domain
figure;
[pxx_manual, f_ps_manual] = pspectrum(Low_level_IIR_Output, fs);
plot(f_ps_manual / 1e3, 10*log10(pxx_manual));
hold on;

% Plot the Un-Filtered signal in frequency domain
[pxx_unfiltered, f_ps] = pspectrum(Signal, fs);
plot(f_ps / 1e3, 10*log10(pxx_unfiltered)); % Convert frequency to kHz and power to dB
xlabel('Frequency (kHz)');
ylabel('Power Spectrum (dB)');
title('Comparison of Filtered and Un-filtered signals in Frequency domain');
legend('Filtered Signal', 'Un-Filtered Signal', 'Location', 'Best');
grid on;

%% Envelope Detection

% Full Wave Rectification of Filtered Signal
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

% Construction of an Ideal Impulse Response
% Specifications
passband_edge = 1e3; % Passband edge frequency in Hz
stopband_edge = 2e3; % Stopband edge frequency in Hz

% Calculate Cutoff Frequency
f_cut = (passband_edge + stopband_edge) / 2; % Midpoint between passband and stopband edges
Fc = f_cut / fs; % Normalized cutoff frequency with respect to sampling frequency

% Calculate Number of Taps (N) using Hamming Window Parameters
delta_f = stopband_edge - passband_edge; % Transition width in Hz
delta_f_normalized = delta_f / fs; % Normalized transition width
N_taps = round(3.3 / delta_f_normalized); % Calculate number of taps
if mod(N_taps, 2) == 0
    N_taps = N_taps + 1; % Ensure the filter length is odd
end

% Generate Ideal Impulse Response
M = (N_taps - 1) / 2;
h_positive = zeros(1, M); % Preallocate for positive indices
for n = 1:M
    h_positive(n) = 2 * Fc * sin(n * 2 * pi * Fc) / (n * 2 * pi * Fc); % Ideal impulse response for LP filter
end
h_ideal = [fliplr(h_positive), 2 * Fc, h_positive]; % Construct full impulse response

% Apply Hamming Window
K_taps = (0:1:N_taps-1)';
window_1stage = 0.54 - (0.46 * (cos(2 * pi * (K_taps /N_taps))));
h_fir = h_ideal .* window_1stage;

% Low-Level FIR Convolution Implementation
tic; % Start timing
y_fir_manual = zeros(size(rectified_signal)); % Initialize output
for q = 1:length(rectified_signal)
    % Convolution sum for FIR filter
    for k = 1:N_taps
        if (q - k + 1) > 0
            y_fir_manual(q) = y_fir_manual(q) + h_fir(k) * rectified_signal(q - k + 1);
        end
    end
end
elapsed_time_1_stage = toc; % Stop timing

% Print the Number of Taps
fprintf('Number of Taps for Single-Stage FIR Filter: %d\n', N_taps);

% Visualize Frequency Response
[H_1stage, W_1stage] = freqz(y_fir_manual, 1, 4096, fs);
figure;
plot(W_1stage / 1e3, 20*log10(abs(H_1stage)));
xlabel('Frequency (kHz)');
ylabel('Magnitude');
title('Frequency Response of Single-Stage FIR Low-Pass Filter');
grid on;

% Plot Filtered Signal in Time Domain
figure;
plot(t_rect, y_fir_manual);
xlabel('Time (s)');
ylabel('Amplitude');
title('Manually Filtered in Time Domain');
grid on;

% Plot Filtered Signal in Frequency Domain
figure;
[pxx_manual_filtered, f_ps_manual_filtered] = pspectrum(y_fir_manual, fs);
plot(f_ps_manual_filtered / 1e3, 10*log10(pxx_manual_filtered));
xlabel('Frequency (kHz)');
ylabel('Power Spectrum (dB)');
title('Filtered Signal Spectrum (Single-Stage FIR)');
grid on;

% Display Timing for Single-Stage Filter
fprintf('Elapsed Time for Single-Stage FIR Filtering (Single stage): %.6f seconds\n', elapsed_time_1_stage);

% Calculate Time Delay and Range of Target
% Convert the time threshold to samples to ignore the initial sonar pulse
ignore_samples = 2100; % Number of samples to ignore outgoing pulse

% Find the peak of the filtered rectified signal after the initial pulse
[~, max_index] = max(y_fir_manual(ignore_samples:end)); % Find peak after the initial 1 ms

% Calculate the Time Delay of the Echo Signal
% Time Delay = Time of Echo - Time of Initial Pulse (ignoring 1 ms)
time_delay = (max_index / fs); % Time delay in seconds

% Calculate Range to Target
speed_of_sound = 1500; % Speed of sound in water (m/s)
range = (time_delay * speed_of_sound) / 2; % Divide by 2 for two-way travel

% Display Results
fprintf('Time Delay: %.6f seconds\n', time_delay);
fprintf('Range to Target: %.2f meters\n', range);

%% Multirate Low Pass Filter

% Decimation factor
decimation_factor = 50;
fs_new = fs / decimation_factor; % New sampling frequency
f_nyquist_new = fs_new / 2; % Nyquist frequency after decimation

% First-Stage Anti-Aliasing Low-Pass Filter
% Specifications
passband_edge_anti = 1.5e3; % Passband edge (1.5 kHz)
stopband_edge_anti = 5e3; % Stopband edge (5 kHz)
transition_width_anti = stopband_edge_anti - passband_edge_anti; % Transition width
f_cutoff_anti = (passband_edge_anti + stopband_edge_anti) / 2; % Cutoff frequency (3.25 kHz)

% Normalized cutoff frequency and transition width
Fc_normalized_anti = f_cutoff_anti / (fs);
delta_f_normalized_anti = transition_width_anti / fs;

% Calculate the number of taps
N_anti = floor(3.3 / delta_f_normalized_anti);
if mod(N_anti, 2) == 0
    N_anti = N_anti + 1; % Ensure odd number of taps
end

% Generate Ideal Impulse Response
M_anti = (N_anti - 1) / 2;
h_positive_anti = zeros(1, M_anti);
for q = 1:M_anti
    h_positive_anti(q) = 2 * Fc_normalized_anti * sin(q * 2 * pi * Fc_normalized_anti) / (q * 2 * pi * Fc_normalized_anti);
end
h_ideal_anti = [fliplr(h_positive_anti), 2 * Fc_normalized_anti, h_positive_anti];

% Apply Hamming Window
K_anti = (0:1:N_anti-1)';
window_anti = 0.54 - (0.46 * (cos(2 * pi * (K_anti /N_anti))));
h_anti = h_ideal_anti .* window_anti;

% Print the Number of Taps
fprintf('Number of Taps for Anti-aliasing Filter: %d\n', N_anti);

% Visualize Frequency Response of First-Stage Filter
[H_anti, W_anti] = freqz(h_anti, 1, 4096, fs);
figure;
plot(W_anti / 1e3, 20*log10(abs(H_anti)));
xlabel('Frequency (kHz)');
ylabel('Magnitude');
title('Frequency Response of First-Stage Anti-Aliasing Filter');
grid on;

% Second-Stage Low-Pass Filter Design
% Specifications
passband_edge = 1e3; % Passband edge frequency in Hz
stopband_edge = 2e3; % Stopband edge frequency in Hz
transition_width = stopband_edge - passband_edge; % Transition width
f_cut = (passband_edge + stopband_edge) / 2;

% Normalized cutoff frequency and transition width
Fc_normalized = f_cut / (fs_new);
transition_width_normalized = transition_width / fs_new;

% Calculate the number of taps
N_stage2 = ceil(3.3 / transition_width_normalized);
if mod(N_stage2, 2) == 0
    N_stage2 = N_stage2 + 1; % Ensure odd number of taps
end

% Generate Ideal Impulse Response
M_stage2 = (N_stage2 - 1) / 2;
h_positive_stage2 = zeros(1, M_stage2);
for r = 1:M_stage2
    h_positive_stage2(r) = 2 * Fc_normalized * sin(r * 2 * pi * Fc_normalized) / (r * 2 * pi * Fc_normalized);
end
h_ideal_stage2 = [fliplr(h_positive_stage2), 2 * Fc_normalized, h_positive_stage2];

% Apply Hamming Window
K_stage2 = (0:1:N_stage2-1)';
window_stage2 = 0.54 - (0.46 * (cos(2 * pi * (K_stage2 /N_stage2))));
h_stage2 = h_ideal_stage2 .* window_stage2;

% Print the Number of Taps
fprintf('Number of Taps for Second-Stage FIR Filter: %d\n', N_stage2);

% Visualize Frequency Response of Second-Stage Filter
[H_stage2, W_stage2] = freqz(h_stage2, 1, 4096, fs_new);
figure;
plot(W_stage2 / 1e3, 20*log10(abs(H_stage2)));
xlabel('Frequency (kHz)');
ylabel('Magnitude');
title('Frequency Response of Second-Stage Anti-Aliasing Filter');
grid on;

% Combined Timing for Both Stages
tic; % Start timing

% Optimized Low-Level Implementation for First-Stage Decimation
first_stage_output = zeros(1, ceil(length(rectified_signal) / decimation_factor));
output_index = 1; % Index for storing decimated output

for n = 1:decimation_factor:length(rectified_signal) % Step by decimation factor
    % Compute convolution for the current decimated sample
    y_temp = 0;
    for k = 1:N_anti
        if (n - k + 1) > 0 % Ensure valid index
            y_temp = y_temp + h_anti(k) * rectified_signal(n - k + 1);
        end
    end
    
    % Store the result in the decimated output
    first_stage_output(output_index) = y_temp;
    output_index = output_index + 1;
end

% Uncomment the line below to use the built-in `upfirdn()` function
% first_stage_output = upfirdn(rectified_signal, h_fir_anti, 1, decimation_factor);

% Low-Level Implementation for Second-Stage Filtering
final_output = zeros(1, length(first_stage_output));
for n = 1:length(first_stage_output)
    y_temp = 0;
    for k = 1:N_stage2
        if (n - k + 1) > 0
            y_temp = y_temp + h_stage2(k) * first_stage_output(n - k + 1);
        end
    end
    final_output(n) = y_temp;
end

% Uncomment the line below to use the built-in `filter()` function
% final_output = filter(h_fir_stage2, 1, first_stage_output);

elapsed_time_multi_rate = toc; % Stop timing

% Plot Spectrum After First Stage
figure;
[pxx_stage1, f_ps_stage1] = pspectrum(first_stage_output, fs_new);
plot(f_ps_stage1 / 1e3, 10 * log10(pxx_stage1));
xlabel('Frequency (kHz)');
ylabel('Power Spectrum (dB)');
title('Frequency Spectrum After First-Stage Decimation');
grid on;

% Plot Spectrum After Second Stage
figure;
[pxx_stage2, f_ps_stage2] = pspectrum(final_output, fs_new);
plot(f_ps_stage2 / 1e3, 10 * log10(pxx_stage2));
xlabel('Frequency (kHz)');
ylabel('Power Spectrum (dB)');
title('Frequency Spectrum After Second-Stage Filtering');
grid on;

% Display Results
fprintf('Elapsed Time for Both Stages (Multi-rate): %.6f seconds\n', elapsed_time_multi_rate);
