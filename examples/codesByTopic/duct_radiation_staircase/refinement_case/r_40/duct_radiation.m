% duct radiation
clear('all'); clc;
%close all;
% chirp
pressures = textread('history_pressures_6r.dat');
signal_in = textread('signal_in.dat');

%pressures = [pressures(:) - mean(pressures)];
pressures = pressures(100:end);
particle_velocity = textread('history_velocities_6r.dat');
particle_velocity = [particle_velocity(:)];
particle_velocity = particle_velocity(100:end);
% a_1 eh o maior / a_2 eh o menor
a_1 = 40;
a_2 = a_1;
cs = 1/sqrt(3);
L = 6*a_1;

particle_velocity = particle_velocity(1:length(pressures));

fft_pressure = fft(pressures);
fft_particle_velocity = fft(particle_velocity);
ZL = fft_pressure./fft_particle_velocity;
frequencies = linspace(0, 1, length(ZL));
frequencies = frequencies';
ka = (2*pi*frequencies*a_1)/cs;
k = ka/a_1;
Zo = 1*cs;%/(pi*a^2);
Zr = Zo*1i*tan(atan(ZL./(1i*Zo))-(ka/a_1)*L);

Rr=(Zr-Zo)./(Zr+Zo);
RL=(ZL-Zo)./(ZL+Zo);
l = (log(-Rr./(abs(Rr))))./(-2*i*k);
%l = (log(-RL./(abs(RL))))./(-2*i*k);
l_line = real(l);
la = l_line/(a_2);
% Refletion Coeff. 
%figure;
%load munt_R.mat
%plot(munt005(:,1),munt005(:,2),'black')
%hold on
%figure;
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/abs_r.fig',1)
hold on;
%plot(ka,abs(RL),'--black');
plot(ka,abs(Rr),'--black');
axis([0 2.5 0 1.5]);
ylabel('Coeficiente de Reflex\E3o, Rr','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
%legend('Analitico','LBM');
hold off

% %% end correction
%figure;
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/loa.fig',1)
% %load munt_loa.mat
% %plot(loa_005(:,1),loa_005(:,2),'black')
hold on
plot(ka,la)
ylabel('End correction, l/a','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
axis([0 1.5 -1 1]);
legend('Analitico','LBM');
% hold off

uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/levine.fig',1)
hold on
plot(ka, real(Zr),'--'); hold on; plot(ka, imag(Zr), 'r--');
%plot(ka, real(ZL),'--'); hold on; plot(ka, imag(ZL), 'r--');
%axis([0 1.82 0 max(abs(ZL))]);
ylabel('Impedance','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
legend('Impedance real part','Impedance imaginary part');
xlim([0 3])

%cd ~/palabos_acoustic/examples/codesByTopic/duct_radiation/
%system('rm tmp/*');
%system('make && time mpirun -n 6 duct_radiation');
