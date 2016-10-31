% get_duct_radiation_parameters: function description
function parameters = get_duct_radiation_parameters(history_pressures_file, history_velocities_file, radius, position)

	pressures = textread(history_pressures_file);
	pressures = [pressures(:)];

	particle_velocity = textread(history_velocities_file);
	particle_velocity = [particle_velocity(:)];

	% filtering signal
	window = hanning(length(pressures));
	pressures = pressures.*window;
	particle_velocity = particle_velocity.*window;

	a_1 = radius;
	cs = 1/sqrt(3);
	L = position;

	a_2 = a_1;
	particle_velocity = particle_velocity(1:length(pressures));
	fft_pressure = fft(pressures);
	fft_particle_velocity = fft(particle_velocity);
	ZL = fft_pressure./fft_particle_velocity;
	frequencies = linspace(0, 1, length(ZL));
	frequencies = frequencies';
	ka = (2*pi*frequencies*a_1)/cs;
	k = ka/a_1;
	Zo = 1*cs;

	Zr = Zo*1i*tan(atan(ZL./(1i*Zo))-(ka/a_1)*L);
	Rr=(Zr-Zo)./(Zr+Zo);

	l = (log(-Rr./(abs(Rr))))./(-2*i*k);
	l_line = real(l);
	la = l_line/(a_2);
	
	parameters{1} = ka;
	parameters{2} = Zr;
	parameters{3} = Rr;
	parameters{4} = la;
