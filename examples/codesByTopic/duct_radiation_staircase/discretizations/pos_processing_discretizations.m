% pos_processing

clear('all'); clc; close all;
cs = 1/sqrt(3);
ka_max = 1.8;
pressures_20 = textread('radius_20/history_pressures_boca.dat');
pressures_20 = [pressures_20(:)];
particle_velocity_20 = textread('radius_20/history_velocities_boca.dat');
particle_velocity_20 = [particle_velocity_20(:)];
a_20 = 20;
fc_20 = (cs*ka_max)/(a_20*2*pi);
B_20 = fir1(60, fc_20, 'low');
pressures_20 = filter(B_20, 1, pressures_20);
particle_velocity_20 = filter(B_20, 1, particle_velocity_20);

pressures_22 = textread('radius_22/history_pressures_boca.dat');
pressures_22 = [pressures_22(:)];
particle_velocity_22 = textread('radius_22/history_velocities_boca.dat');
particle_velocity_22 = [particle_velocity_22(:)];
a_22 = 22;

pressures_24 = textread('radius_24/history_pressures_boca.dat');
pressures_24 = [pressures_24(:)];
particle_velocity_24 = textread('radius_24/history_velocities_boca.dat');
particle_velocity_24 = [particle_velocity_24(:)];
a_24 = 24;

pressures_40 = textread('radius_40/history_pressures_boca.dat');
pressures_40 = [pressures_40(:)];
particle_velocity_40 = textread('radius_40/history_velocities_boca.dat');
particle_velocity_40 = [particle_velocity_40(:)];
a_40 = 40;
fc_40 = (cs*ka_max)/(a_40*2*pi);
B_40 = fir1(60, fc_40, 'low');
pressures_40_f = filter(B_40, 1, pressures_40);
particle_velocity_40_f = filter(B_40, 1, particle_velocity_40);

fft_pressure_20 = fft(pressures_20);
fft_particle_velocity_20 = fft(particle_velocity_20);
ZL_20 = fft_pressure_20./fft_particle_velocity_20;

fft_pressure_22 = fft(pressures_22);
fft_particle_velocity_22 = fft(particle_velocity_22);
ZL_22 = fft_pressure_22./fft_particle_velocity_22;

fft_pressure_24 = fft(pressures_24);
fft_particle_velocity_24 = fft(particle_velocity_24);
ZL_24 = fft_pressure_24./fft_particle_velocity_24;


fft_pressure_40 = fft(pressures_40);
fft_pressure_40_f = fft(pressures_40_f);

fft_particle_velocity_40 = fft(particle_velocity_40);
fft_particle_velocity_40_f = fft(particle_velocity_40_f);

ZL_40 = fft_pressure_40./fft_particle_velocity_40;
ZL_40_f = fft_pressure_40_f./fft_particle_velocity_40_f;

frequencies = linspace(0, 1, length(ZL_20));
frequencies = frequencies';
ka_20 = (2*pi*frequencies*a_20)/cs;
ka_24 = (2*pi*frequencies*a_24)/cs;
ka_22 = (2*pi*frequencies*a_22)/cs;
ka_40 = (2*pi*frequencies*a_40)/cs;
k = ka_22/a_22;
Zo = 1*cs;%/(pi*a^2);

%A_20=pi*a_20^2;
%A_22=pi*a_22^2;
%A_24=pi*a_24^2;

uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/levine.fig',1)
hold on
%figure;
%plot(ka, real(Zr),'--'); hold on; plot(ka, imag(Zr), 'r--');
plot(ka_40, real(ZL_40/Zo),'--', 'Color',[0.5,0,0.5]); hold on; plot(ka_40, imag(ZL_40/Zo), '--', 'Color',[0.5,0,0.5]);
%plot(ka_22, real(ZL_22/Zo),':', 'Color',[0.5,0.5,0]); hold on; plot(ka_22, imag(ZL_22/Zo), '--', 'Color',[0.5,0.5,0]);
%plot(ka_24, real(ZL_24/Zo),':', 'Color',[0,0.5,0.5]); hold on; plot(ka_24, imag(ZL_24/Zo), '--', 'Color',[0,0.5,0.5]);
plot(ka_40, real(ZL_40_f/Zo),'--'); hold on; plot(ka_40, imag(ZL_40_f/Zo), 'r--');
%axis([0 1.82 0 max(abs(ZL))]);
ylabel('Impedance','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
legend('Impedance real Levine','Impedance imaginary Levine', ... 
	'Impedance real part-10','Impedance imaginary part-10', ...
	'Impedance real part-20','Impedance imaginary part-20');
%xlim([0 2])
axis([0 3 0 1.2])