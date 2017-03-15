% multiple_microphone_method
close all; clear('all'); clc;

%% visc = (-1/2 + 1/omega)*cs2;
%%  Re = ((M*cs)*40)/visc

%% WAVE DECOMPOSITION
a = 20;
cs = 1/sqrt(3);
M = -0.15;
ka_max_excitation = 2.5;
% estou_esperancoso_no_Mats_Abom_Munt
% testando_um_dominio_menor
% o_primeiro_cara_sugado
signal_in = textread('o_primeiro_cara_sugado/signal_in.dat');
signal_in = signal_in - rms(signal_in);
frequencies = linspace(0,1,length(signal_in));
ka = (2*pi*frequencies*a)/cs;
[picos, locais] = findpeaks(abs(fft(signal_in)));
ka_picos = ka(locais);
ka_picos = ka_picos(find(ka_picos <= ka_max_excitation));

%iT <= maxT_final_source /*&& iT > maxT/2)
pressures = textread('o_primeiro_cara_sugado/system_abom_measurement_data/history_pressures_system_abom_measurement_data.dat');
length_duct = 3*(6*a);
nz = length_duct + 6*a + 30;
maxT = 2*(2^13 + nz*sqrt(3));
%maxT_final_source = maxT - nz*sqrt(3);
pressures = pressures(round(maxT/2):end, :);

% Setting up mics positions
[size_1 size_2] = size(pressures);
%final_mic_position = size_2 - 13;
final_mic_position = size_2 - 13 + a/4;
start_mic_position = final_mic_position - round(a*4);
%start_mic_position = 100;
positions = start_mic_position:final_mic_position;
reference_position_mic = final_mic_position;

pressures_frequencies(1:length(ka_picos), 1:length(positions)) = 0;
for position = 1:length(positions)
	signal = pressures(:, positions(position));
	window = hanning(length(signal));
	signal = window.*signal;
	ka_signal = (2*pi*linspace(0,1,length(signal))*a)/cs;
	fft_signal = fft(signal);
	
	for pico = 1:length(ka_picos)
		[value slot] = min(abs(ka_signal - ka_picos(pico)));
		pressures_frequencies(pico, position) = fft_signal(slot);
	end

end

nFreq = length(ka_picos);
for f = 1:nFreq
    ka = ka_picos(f); % 500 to 3000 Hz (20 step)
    frequency = (ka*cs)/(a*2*pi);
    w = frequency*2*pi;
    
    k1 = w/cs;
    ki(f) = k1/(1 + M);
    kr(f) = k1/(1 - M);
    
    pressure2 = pressures_frequencies(f, :).'; %pressure on all mics after test section (complex numbers)
    [pplus2,pminus2] = waveDecomposition(pressure2,positions,reference_position_mic,ki(f),kr(f));
    Re(f) = pminus2/pplus2;
end

%/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/abs_r.fig
%panhuis_code/abs_r_015.fig
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/abs_r.fig',1)
hold on;
plot(ka_picos,smooth(abs(Re)),'*black');
axis([0 2.5 0 1.2]);

k_picos = ka_picos/a;
l = (log(-Re./(abs(Re))))./(-2*i*k_picos);
l_line = real(l);
la = l_line/(a);

uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/loa.fig',2)
hold on;
plot(ka_picos,smooth(la), '*black')
ylabel('End correction, l/a','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
axis([0 1.5 0 0.8]);

saveas(1, 'escoamento_sugado_abs_r', 'png');
saveas(2, 'escoamento_sugado_loa', 'png');