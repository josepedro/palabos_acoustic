% We can do that! Corollary Howe

clear('all'); close all; clc;

total_time = 482;
initial_time = 0;
times_period = initial_time:total_time-1;
radius = 20;
diameter = 2*radius;
length_duct = 3*(3*diameter);
distance_from_start = 0; 
Nz = 1 + length_duct + diameter - distance_from_start;

file_name = 'analisando_howe/max_ka/';

sentinella = textread(strcat(file_name,'howe_corollary_0_upright.dat'));
sentinella = sentinella(1:end-1);
sentinella = reshape(sentinella, Nz, length(sentinella)/Nz);
sizes = size(sentinella);

% abrindo arquivos e colocando nas matrizes
velocities_axial(1:sizes(1), 1:sizes(2), 1:length(times_period)) = 0; 
velocities_upright(1:sizes(1), 1:sizes(2), 1:length(times_period)) = 0;
veltorial_product_axial = velocities_upright;
veltorial_product_upright = velocities_upright; 
mean_velocities_axial(1:sizes(1), 1:sizes(2)) = 0;
mean_velocities_upright = mean_velocities_axial;
for time = 1:length(times_period)
	sentinella = textread(strcat(file_name, 'howe_corollary_', num2str(time - 1), '_axial.dat'));
	sentinella = sentinella(1:end-1);
	sentinella = reshape(sentinella, Nz, length(sentinella)/Nz);
	velocities_axial(:,:,time) = sentinella;
	mean_velocities_axial = mean_velocities_axial + velocities_axial(:,:,time); 

	sentinella = textread(strcat(file_name, 'howe_corollary_', num2str(time - 1), '_upright.dat'));
	sentinella = sentinella(1:end-1);
	sentinella = reshape(sentinella, Nz, length(sentinella)/Nz);
	velocities_upright(:,:,time) = sentinella;
	mean_velocities_upright = mean_velocities_upright + velocities_upright(:,:,time);

	% calculando o rotacional
	omega = curl(velocities_axial(:,:,time), velocities_upright(:,:,time));

	% calculando produto vetorial
	for z = 1:sizes(1)
		for x = 1:sizes(2)
			alphA = omega(z, x); 
			betA = velocities_axial(z,x,time);
			gammA = velocities_upright(z,x,time);
			 
			cross_product = cross([0,0,alphA],[betA,gammA,0]);

			veltorial_product_axial(z,x,time) = cross_product(1);
			veltorial_product_upright(z,x,time) = cross_product(2);  
		end		
	end
end

% calculando as medias
mean_velocities_upright = mean_velocities_upright/length(times_period);
mean_velocities_axial = mean_velocities_axial/length(times_period);

% calulando o resultado final
acoustic_energy(1:sizes(1), 1:sizes(2), 1:length(times_period)) = 0;
for time = 1:length(times_period)
	for z = 1:sizes(1)
		for x = 1:sizes(2)
			alphA = veltorial_product_axial(z,x,time);
			betA = veltorial_product_upright(z,x,time);
			Uac_axial = velocities_axial(z,x,time) - mean_velocities_axial(z,x);  
			Uac_upright = velocities_upright(z,x,time) - mean_velocities_upright(z,x); 			 
			acoustic_energy(z,x,time) = dot([alphA betA],[Uac_axial Uac_upright]);
		end		
	end
end
