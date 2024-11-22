%% Task 1
% Time and frequency domain analysis

% Load the Signal
[sonar_signal, fs] = audioread('Shreyas Ravi.wav');

% Inspect the outgoing Sonar Pulse and rest of the recieved signal
t = (0:length(sonar_signal)-1) / fs; 
Fs = 500000;
Ts = 1/Fs;
Tmax = 1;
t = 0:Ts:Tmax-Ts;      % Generate discrete time values (nTs)

figure; 
plot(t, Sonar_Signal(1:500));    % Plot
xlabel('Time (Sec)'); ylabel('Amplitude'); title('Sonar Signal');
grid on;