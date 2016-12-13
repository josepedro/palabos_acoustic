%% get_anechoic_parameters: function description
function parameters = get_duct_radiation_parameters_boca(pressures, history_velocities_x, history_velocities_y, history_velocities_z, size_buffer)

	pressures = [pressures(:)];
	%size(pressures)
	%figure(1); plot(pressures);


	particle_velocity = history_velocities_z;
	particle_velocity = [particle_velocity(:)];
	particle_velocity = particle_velocity(1:length(pressures));
	%size(particle_velocity)
	%figure(2); plot(particle_velocity);

	% filtering signal
	window = hanning(length(pressures));
	pressures = pressures.*window;
	particle_velocity = particle_velocity.*window;

	a_1 = size_buffer;
	cs = 1/sqrt(3);

	a_2 = a_1;
	particle_velocity = particle_velocity(1:length(pressures));
	fft_pressure = fft(pressures);
	fft_particle_velocity = fft(particle_velocity);
	ZL = fft_pressure./fft_particle_velocity;
	frequencies = linspace(0, 1, length(ZL));
	frequencies = frequencies';
	ka = (2*pi*frequencies*a_1)/cs;
	k = ka/a_1;
	Zo = 1*cs;%/(pi*a^2);
	RL=(ZL-Zo)./(ZL+Zo);

	l = (log(-RL./(abs(RL))))./(-2*i*k);
	l_line = real(l);
	la = l_line/(a_2);

	
	parameters{1} = ka;
	parameters{2} = ZL;
	parameters{3} = RL;
	parameters{4} = la;
