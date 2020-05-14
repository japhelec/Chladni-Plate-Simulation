clear all;
filename = '2025_2Dfreq.txt';
ts = 0.001;
s = textread(filename);
fs = 1 / ts;
f = linspace(-fs/2, fs/2, length(s));
S = fftshift(fft(s));
plot(f, abs(S));
