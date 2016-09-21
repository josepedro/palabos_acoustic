% duct radiation

% chirp
pressures = textread('tmp/history_pressures_3r.dat');
pressures = [pressures(:)];
particle_velocity = textread('tmp/history_velocities_3r.dat');
particle_velocity = [particle_velocity(:)];
a = 6;
cs = 1/sqrt(3);
L = a*3; 

particle_velocity = particle_velocity(1:length(pressures));

fft_pressure = fft(pressures);
fft_particle_velocity = fft(particle_velocity);
ZL = fft_pressure./fft_particle_velocity;
%ka_max = (2*pi*a)/cs;
frequencies = linspace(0, 1, length(ZL));
frequencies = frequencies';
ka = (2*pi*frequencies*a)/cs;
k=ka/a;
%Zo = rms(dens_filter)*cs;
Zo = 1*cs;
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

figure;
plot(ka, real(Zr)); hold on; plot(ka, imag(Zr), 'r');
ylabel('Impedance','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
legend('Impedance real part','Impedance imaginary part');


%cd ~/palabos_acoustic/examples/codesByTopic/duct_radiation/
%system('rm tmp/*');
%system('make && time mpirun -n 6 duct_radiation');
