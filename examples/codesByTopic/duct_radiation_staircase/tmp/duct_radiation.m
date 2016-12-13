% duct radiation
clear('all'); clc;
%close all;
% chirp
pressures = textread('history_pressures_boca.dat');
pressures = [pressures(:) - mean(pressures)];
%pressures = pressures(5000:end);
particle_velocity = textread('history_velocities_boca.dat');
particle_velocity = [particle_velocity(:)];
%particle_velocity = particle_velocity(5000:end);
a = 18;
cs = 1/sqrt(3);
L = 26;

particle_velocity = particle_velocity(1:length(pressures));

fft_pressure = fft(pressures);
fft_particle_velocity = fft(particle_velocity);
ZL = fft_pressure./fft_particle_velocity;
%ZL = ZL;
%ka_max = (2*pi*a)/cs;
frequencies = linspace(0, 1, length(ZL));
frequencies = frequencies';
ka = (2*pi*frequencies*a)/cs;
k=ka/a;
%Zo = rms(dens_filter)*cs;
Zo = 1*cs;%/(pi*a^2);
Zr = Zo*1i*tan(atan(ZL./(1i*Zo))-(ka/a)*L);


Rr=(Zr-Zo)./(Zr+Zo);
RL=(ZL-Zo)./(ZL+Zo);
l=(1./k).*atan(Zr/(1i*Zo));
la=real(l/a);
% Refletion Coeff. 
%figure;
%load munt_R.mat
%plot(munt005(:,1),munt005(:,2),'black')
%hold on
%figure;
% plot(ka,abs(Rr),'--black')
% axis([0 1.5 0 5]);
% ylabel('Coeficiente de Reflex\E3o, Rr','FontSize',20);
% xlabel('Numero de Helmholtz, ka','FontSize',20);
% %legend('Analitico','LBM');
% hold off

% %% end correction
% figure;
% %load munt_loa.mat
% %plot(loa_005(:,1),loa_005(:,2),'black')
% %hold on
% plot(ka,la)
% ylabel('End correction, l/a','FontSize',20);
% xlabel('Numero de Helmholtz, ka','FontSize',20);
% axis([0 1.5 -1 1]);
% %legend('Analitico','LBM');
% hold off
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/levine.fig',1)
hold on

%figure;
%plot(ka, real(Zr),'--'); hold on; plot(ka, imag(Zr), 'r--');
plot(ka, real(ZL/Zo),'--'); hold on; plot(ka, imag(ZL/Zo), 'r--');
%axis([0 1.82 0 max(abs(ZL))]);
ylabel('Impedance','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
legend('Impedance real part','Impedance imaginary part');
xlim([0 2])


%cd ~/palabos_acoustic/examples/codesByTopic/duct_radiation/
%system('rm tmp/*');
%system('make && time mpirun -n 6 duct_radiation');
