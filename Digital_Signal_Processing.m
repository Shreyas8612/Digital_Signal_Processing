% Time and frequency domain analysis

% Load the Signal
[sonar_signal, fs] = audioread('Sonar_Signal.wav');

% Inspect the outgoing Sonar Pulse and rest of the recieved signal
t = (0:length(sonar_signal)-1) / fs; % Time Vector
sonar_signal_samples = round(1e-3 * fs); % Number of samples for the first 1ms

% Plot the outgoing pulse
figure;
plot(t(1:sonar_signal_samples), sonar_signal(1:sonar_signal_samples));
xlabel('Time (s)');
ylabel('Amplitude');
title('Outgoing SOnar Pulse in Time Domain');
grid on;

% Plot the entire received signal
figure;
plot(t, sonar_signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Received Sonar Signal in Time Domain');
grid on;

% Frequency Domain Analysis using fft() fast fourier transforms
N = length(sonar_signal);
fft_signal = fft(sonar_signal); % Compute FFT
f = (0:N-1) * (fs / N); % Frequency vector for FFT

% Normalization and logarithmic scaling
fft_magnitude = abs(fft_signal) / N;
fft_magnitude_db = 20 * log10(fft_magnitude); % Converty to dB

% Plot +ve half of the spectrum
figure;
plot(f(1:N/2), fft_magnitude_db(1:N/2));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Domain Representation of Sonar Signal');
grid on;

% Identify the carrier frequency
[~, peak_index] = max(fft_magnitude(1:N/2));
carrier_frequency = f(peak_index);
fprintf('The carrier frequency is %.1f kHz\n', carrier_frequency / 1e3);