% find frequency resonator
close all
clear('all');
pressures = textread('velocities_y.dat');
pressures = pressures - mean(pressures);
window = hanning(length(pressures));
pressures = pressures.*window;
velocities = textread('velocities_y.dat');
velocities = velocities - mean(velocities);
velocities = velocities.*window;
%impedance = fft(pressures)./fft(velocities);


N=length(pressures);
windowLenght=N/80;
%Sxx = cpsd(pressures,pressures,hanning(windowLenght),N/90,N,'twosided');
%Sxx = cpsd(pressures,pressures);
frequencies = linspace(0, 1, length(pressures));
figure(1);
%plot(frequencies,abs(Sxx));
plot(frequencies, abs(fft(pressures)));
xlim([0 0.1]);
%ylim([0 0.08])
