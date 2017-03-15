% test function to filter kas
clear all; close all; clc;

a = 20;
cs = 1/sqrt(3);
signal_in = textread('signal_in.dat');
frequencies = linspace(0,1,length(signal_in));
ka = (2*pi*frequencies*a)/cs;
[picos, locais] = findpeaks(abs(fft(signal_in)));
ka_picos = ka(locais);

% escolha seu pico
it = [0:5000];
phase = (ka_picos(1)*cs*it)/a;
sin_1 = sin(phase);
phase = (ka_picos(2)*cs*it)/a;
sin_2 = sin(phase);

% convoluindo os sinais
result_1 = conv(signal_in, sin_1);
result_2 = conv(signal_in, sin_2);

result_total = result_1 + result_2;

figure; plot(abs(fft(signal_in)))

figure; plot(abs(fft(result_total)))


