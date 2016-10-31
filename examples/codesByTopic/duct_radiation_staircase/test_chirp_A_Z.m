% test chirp A_Z
%clear all; close all; clc;

time_total = 2^13;
chirp_signal(1:time_total) = 0;
amplitude = 1;
cs = 1/sqrt(3);
a = 20;
frequencies = linspace(0, 1, time_total);
ka = (2*pi*frequencies*a)/cs;
%n_sin = 1;
ka_max = 2.5;
total_n = 20;
for time = 1:time_total
	for n_sin = 1:total_n
		%random_phase = rand(1,1)*pi;
		random_phase = 0;
		chirp_signal(time) = chirp_signal(time) + amplitude*sin(random_phase+(n_sin*(ka_max/total_n)*cs*time)/(a));
	end
end
chirp_signal = chirp_signal/max(abs(chirp_signal));

figure(1);
plot(chirp_signal);
%hold on;
%plot(signal_in - mean(signal_in), 'r'); hold off;

figure(2);
plot(ka, abs(fft(chirp_signal)));
%hold on;
%fft_signal_in = fft(signal_in - mean(signal_in));
%fft_signal_in(1:2) = 0;
%plot(ka, abs(fft_signal_in), 'r');
xlim([0 2.6]);


% total_time = time_total;
% times = 0 : total_time - 1;
% initial_frequency = 0*cs/(2*pi*a);
% frequency_max_lattice = ka_max*cs/(2*pi*a);
% source_chirp = chirp(times, ... 
% initial_frequency, times(end), frequency_max_lattice);

% figure(3);
% plot(source_chirp);

% figure(4);
% plot(ka, abs(fft(source_chirp)));
% xlim([0 1.9]);