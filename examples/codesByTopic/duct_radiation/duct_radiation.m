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
figure; plot(source_chirp); hold on; plot(chirp_hand);


%cd ~/palabos_acoustic/examples/codesByTopic/duct_radiation/
%system('rm tmp/*');
%system('make && time mpirun -n 6 duct_radiation');
