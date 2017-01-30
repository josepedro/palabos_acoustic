% find frequency resonator
close all
clear('all');
pressures_in = textread('pressures_in.dat');
pressures_in = pressures_in - rms(pressures_in);
pressures_out = textread('pressures_out.dat');
pressures_out = pressures_out - rms(pressures_out);

A = 1;
B = 20;
windowLength = length(pressures_in)/A;
Sxx_in = cpsd(pressures_in, pressures_in, hann(windowLength), windowLength/B, length(pressures_in), 1, 'twosided');
psdx_in = abs(Sxx_in(1:length(pressures_in)/2+1));
fft_pressures_in = sqrt(psdx_in*1*length(pressures_in));

windowLength = length(pressures_out)/A;
Sxx_out = cpsd(pressures_out, pressures_out, hann(windowLength), windowLength/B, length(pressures_out), 1, 'twosided');
psdx_out = abs(Sxx_out(1:length(pressures_out)/2+1));
fft_pressures_out = sqrt(psdx_out*1*length(pressures_out));


window = hanning(length(pressures_in));
pressures_in = pressures_in.*window;
pressures_out = pressures_out.*window;

transmission_coefficient_2 = (abs(fft_pressures_out./fft_pressures_in)).^2;
transmission_loss_2 = 10*log10(1./transmission_coefficient_2);

transmission_coefficient_1 = (abs(fft(pressures_out)./fft(pressures_in))).^2;
transmission_loss_1 = 10*log10(1./transmission_coefficient_1);

frequencies_2 = linspace(0, 1, length(transmission_coefficient_2));
frequencies_1 = linspace(0, 1, length(transmission_coefficient_1));
figure(1);
%plot(frequencies,abs(Sxx));smooth(real(Zr), 7)
plot(frequencies_1, transmission_loss_1, 'r');
hold on;
plot(frequencies_1, smooth(transmission_loss_1, 30));
xlim([0 0.5]);
%ylim([0 0.08])
