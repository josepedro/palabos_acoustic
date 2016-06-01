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

compare_pressure_theta_10 = open('pressures_theta.mat');
compare_pressure_theta_10 = compare_pressure_theta_10.pressures_theta;

%figure; polar(theta, pressures_theta, '--');
%hold on; polar(theta, compare_pressure_theta_10, 'red');

%save(filename,variables);

figure; polar(theta, pressures_theta, '--');
hold on; polar(theta, DNS_vortex_shedding(:,1), 'red');
polar(Curle_vortex_shedding(:,2)*(pi/180), Curle_vortex_shedding(:,1), 'black');
title( ... 
'Directivity', ... 
'Interpreter','latex','FontSize',16);
xlabel('Angle [Degree]','Interpreter','latex','FontSize',16); 
ylabel('Normalized Magnitude','Interpreter','latex','FontSize',16);
legend('LBM and FW-H Analogy - Radius Surface = 2','DNS', 'DNS and Curle Analogy');

%title( ... 
%'Time histories of the fluctuation pressure', ... 
%'Interpreter','latex','FontSize',16);
%xlabel('Time step [Lattice unit]','Interpreter','latex','FontSize',16); 
%ylabel('Pressure variation','Interpreter','latex','FontSize',16);

% plot(a.pressure_points_complete_partial);
% axis([0 1001 min(a.pressure_points_complete_partial) max(a.pressure_points_complete_partial)]);
% title( ... 
% 'Time histories of the fluctuation pressure', ... 
% 'Interpreter','latex','FontSize',16);
% xlabel('Time step [Lattice unit]','Interpreter','latex','FontSize',16); 
% ylabel('Pressure variation','Interpreter','latex','FontSize',16);

% frequencies =  linspace(0,1,length(a.pressure_points_complete_partial));
% a.pressure_points_complete_partial = a.pressure_points_complete_partial - mean(a.pressure_points_complete_partial);
% figure; plot(frequencies, abs(fft(a.pressure_points_complete_partial))/max(abs(fft(a.pressure_points_complete_partial)))); 
% axis([0 0.06 0 1.2]);
% title( ... 
% 'Spectrum of Frequencies', ... 
% 'Interpreter','latex','FontSize',16);
% xlabel('Frequencies [Lattice Unit]','Interpreter','latex','FontSize',16); 
% ylabel('Normalized Magnitude','Interpreter','latex','FontSize',16);