% duct_radiation_2

clear('all'); clc;

pressures_1 = textread('history_pressures_3r.dat');
pressures_1 = [pressures_1(:)];
pressures_2 = textread('history_pressures_4r.dat');
pressures_2 = [pressures_2(:)];

particle_velocity_1 = textread('history_velocities_3r.dat');
particle_velocity_1 = [particle_velocity_1(:)];
particle_velocity_2 = textread('history_velocities_4r.dat');
particle_velocity_2 = [particle_velocity_2(:)];

L_1 = 45;
L_2 = 25;
radius = 18;
cs = 1/sqrt(3);
frequencies = linspace(0, 1, length(pressures_1));
frequencies = frequencies';
ka = (2*pi*frequencies*radius)/cs;
k = ka/radius;
rho0 = 1;

B_1 = pressures_1 - pressures_2.*exp(i*k*(L_1 - L_2));
B_2 = exp(i*k*(2*L_1 - L_2)) + exp(i*k*L_2);
B = B_1./B_2;

A_1 = pressures_1 - B.*exp(i*k*L_1);
A_2 = exp(-i*k*L_1);
A = A_1./A_2;

pressures_0 = A + B;
particle_velocity_0 = (A - B)/(cs*rho0);
fft_pressures_0 = fft(pressures_0);
fft_particle_velocity_0 = fft(particle_velocity_0);

Zr = fft_pressures_0./fft_particle_velocity_0;

uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/levine.fig',1)
hold on

plot(ka, real(Zr),'--'); hold on; plot(ka, imag(Zr), 'r--');
