%% Time and frequency domain analysis

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
[pxx, f_ps] = pspectrum(sonar_signal, fs);
plot(f_ps / 1e3, 10*log10(pxx)); % Convert frequency to kHz and power to dB
xlabel('Frequency (kHz)');
ylabel('Power Spectrum (dB)');
title('Frequency Spectrum');
grid on;

% Frequency Domain Analysis using fft() fast fourier transforms
N = length(sonar_signal);
K = (0:N-1)';

% Apply Hanning window
H = 0.5 * (1-cos(2*pi*K/N-1));
Windowed_Signal = sonar_signal .* H;

% Perform FFT on windowed and non-windowed signals
fft_signal_window = fft(Windowed_Signal);
fft_signal = fft(sonar_signal); % Compute FFT
f = (0:N-1) * (fs / N); % Frequency associated with each bin.

% Normalization and logarithmic scaling
fft_magnitude = abs(fft_signal) / N;
fft_magnitude_db = 20 * log10(fft_magnitude); % Convert to dB

% Plot FFT magnitude for original and windowed signals
figure;
plot(f(1:N/2) / 1e3, abs(fft_signal(1:N/2)), 'b');
hold on;
plot(f(1:N/2) / 1e3, abs(fft_signal_window(1:N/2)), 'r');
xlabel('Frequency (kHz)');
ylabel('Magnitude');
title('Frequency Domain Representation of Sonar Signal (Rectangular vs Hanning Window)');
legend('Rectangular Window', 'Hanning Window');
grid on;

% Identify the carrier frequency
[~, peak_index] = max(fft_magnitude(1:N/2));
carrier_frequency = f(peak_index);
fprintf('The carrier frequency is %.1f kHz\n', carrier_frequency / 1e3);