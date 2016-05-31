% script to calculate directivity using FW-H surface
clear('all');
close all;

% getting data from files
sfwh_pressure = open('sfwh_pressure.dat');
sfwh_pressure = sfwh_pressure.sfwh_pressure;
sfwh_velocity_x = open('sfwh_velocity_x.dat');
sfwh_velocity_x = sfwh_velocity_x.sfwh_velocity_x;
sfwh_velocity_y = open('sfwh_velocity_y.dat');
sfwh_velocity_y = sfwh_velocity_y.sfwh_velocity_y;
% prepare matrix ffowcs-williams H. surface
sfwh_matrix(:,:,1) = sfwh_pressure;
sfwh_matrix(:,:,2) = sfwh_velocity_x;
sfwh_matrix(:,:,3) = sfwh_velocity_y;

% open benchmarks data
DNS_vortex_shedding = open('DNS_vortex_shedding.csv');
DNS_vortex_shedding = DNS_vortex_shedding.DNS_vortex_shedding;
DNS_vortex_shedding(:,1) = DNS_vortex_shedding(:,1)/max(DNS_vortex_shedding(:,1)); 
Curle_vortex_shedding = open('Curle_vortex_shedding.csv');
Curle_vortex_shedding = Curle_vortex_shedding.Curle_vortex_shedding;
Curle_vortex_shedding(:,1) = Curle_vortex_shedding(:,1)/max(Curle_vortex_shedding(:,1));

f = 0.006; % frequency peak
%theta = (0:10:360)*pi/180;
theta = DNS_vortex_shedding(:,2)*pi/180; % theta in radians
delta_x = 0.004177e-3;
size_square = 2;
radius = 75*(size_square/delta_x);
%radius = 75*2;
pressures_theta(1:length(theta)) = 0;
U = 0.2/sqrt(3);
V = 0;

for point = 1:length(theta)
	disp((point/length(theta)*100));
	[x y] = pol2cart(theta(point), radius);
	pesc(1) = x;
	pesc(2) = y;
	pressures_theta(point) = fwh(sfwh_matrix,f,U,V,pesc);
end
pressures_theta = abs(pressures_theta)/(U^2);
pressures_theta = pressures_theta/max(pressures_theta);
pressures_theta = pressures_theta';

figure; polar(theta, pressures_theta, '*');
hold on; polar(theta, DNS_vortex_shedding(:,1), 'red');
polar(Curle_vortex_shedding(:,2)*(pi/180), Curle_vortex_shedding(:,1), 'black');
title( ... 
'Directivity', ... 
'Interpreter','latex','FontSize',16);
xlabel('Angle [Degree]','Interpreter','latex','FontSize',16); 
ylabel('Normalized Magnitude','Interpreter','latex','FontSize',16);
legend('Lattice Boltzmann and FW-H Analogy','DNS', 'DNS and Curle Analogy');