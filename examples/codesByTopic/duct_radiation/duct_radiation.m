% duct radiation

% chirp
a = 20;
cs = 1/sqrt(3);
final_time_chirp = 10000;
times = 0 : final_time_chirp - 1;
ka_min = 0;
ka_max = 1.82;
initial_frequency = ka_min*cs/(2*pi*a);
frequency_max_lattice = ka_max*cs/(2*pi*a);
source_chirp = chirp(times, ... 
initial_frequency, times(end), frequency_max_lattice);

variation_frequency = (frequency_max_lattice - initial_frequency)/length(times);
frequency_function = initial_frequency*times + (variation_frequency*times.^2)/2;
phase = 2*pi*frequency_function;
chirp_hand = sin(phase);
%figure; plot(source_chirp); hold on; plot(chirp_hand);

pressures = textread('tmp/history_pressures.dat');
particle_velocity = textread('tmp/history_velocities_x.dat');
fft_pressure = fft(pressures);
fft_particle_velocity = fft(particle_velocity);
ZL = fft_pressure./fft_particle_velocity;
%ka_max = (2*pi*a)/cs;
ka = linspace(0, ka_max, length(fft_pressure));
ka = ka';
k=ka/a;
L = 125; 
%Zo = rms(dens_filter)*cs;
Zo = 1*cs;
Zr = Zo*1i*tan(atan(ZL./(1i*Zo))-(ka/a)*L);
Rr=(Zr-Zo)./(Zr+Zo);
RL=(ZL-Zo)./(ZL+Zo);
l=(1./k).*atan(Zr/(1i*Zo));
la=real(l/a);
% Refletion Coeff. 
figure;
%load munt_R.mat
%plot(munt005(:,1),munt005(:,2),'black')
%hold on
%figure;
plot(ka,abs(Rr),'--black')
axis([0 1.5 0 5]);
ylabel('Coeficiente de Reflex\E3o, Rr','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
%legend('Analitico','LBM');
hold off

%% end correction
figure;
%load munt_loa.mat
%plot(loa_005(:,1),loa_005(:,2),'black')
%hold on
plot(ka,la)
ylabel('End correction, l/a','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
axis([0 1.5 -1 1]);
%legend('Analitico','LBM');
hold off


%cd ~/palabos_acoustic/examples/codesByTopic/duct_radiation/
%system('rm tmp/*');
%system('make && time mpirun -n 6 duct_radiation');
