% duct_radiation_2

clear('all'); clc;

pressures_1 = textread('history_pressures_6r.dat');
pressures_1 = [pressures_1(:)];
pressures_2 = textread('history_pressures_boca.dat');
pressures_2 = [pressures_2(:)];

a = 20;
const = 3;
L_1 = 6*a + const;
L_2 = const;
cs = 1/sqrt(3);
frequencies = linspace(0, 1, length(pressures_1));
frequencies = frequencies';
ka = (2*pi*frequencies*a)/cs;
k = ka/a;
rho0 = 1;

H_12 = (fft(pressures_2))./(fft(pressures_1));
s = L_1 - L_2;
H_I = exp(-i*k*s);
H_R = exp(i*k*s);
A = H_12 - H_I;
B = H_R - H_12;
r = ((A./B)).*exp(2*i*k*L_1);

Zr = (1 + r)./(1 - r);
%Zr = fft_pressures_0./fft_particle_velocity_0;

uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/levine.fig',1)
hold on

plot(ka, real(Zr),'--'); hold on; plot(ka, imag(Zr), 'r--');
xlim([0 1.8])

