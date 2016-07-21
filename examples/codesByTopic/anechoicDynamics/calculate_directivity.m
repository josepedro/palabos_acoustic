% script to calculate directivity using FW-H surface
clear('all');
close all;

% getting data from files
sfwh_pressure_10 = open('sfwh_pressure_10.dat');
sfwh_pressure_10 = sfwh_pressure_10.sfwh_pressure_10;
sfwh_pressure_2 = open('sfwh_pressure_2.dat');
sfwh_pressure_2 = sfwh_pressure_2.sfwh_pressure_2;
sfwh_pressure_50 = open('sfwh_pressure_50.dat');
sfwh_pressure_50 = sfwh_pressure_50.sfwh_pressure_50;
sfwh_pressure_100 = open('sfwh_pressure_100.dat');
sfwh_pressure_100 = sfwh_pressure_100.sfwh_pressure_100;
sfwh_pressure_300 = open('sfwh_pressure_300.dat');
sfwh_pressure_300 = sfwh_pressure_300.sfwh_pressure_300;

sfwh_velocity_x_10 = open('sfwh_velocity_x_10.dat');
sfwh_velocity_x_10 = sfwh_velocity_x_10.sfwh_velocity_x_10;
sfwh_velocity_x_2 = open('sfwh_velocity_x_2.dat');
sfwh_velocity_x_2 = sfwh_velocity_x_2.sfwh_velocity_x_2;
sfwh_velocity_x_50 = open('sfwh_velocity_x_50.dat');
sfwh_velocity_x_50 = sfwh_velocity_x_50.sfwh_velocity_x_50;
sfwh_velocity_x_100 = open('sfwh_velocity_x_100.dat');
sfwh_velocity_x_100 = sfwh_velocity_x_100.sfwh_velocity_x_100;
sfwh_velocity_x_300 = open('sfwh_velocity_x_300.dat');
sfwh_velocity_x_300 = sfwh_velocity_x_300.sfwh_velocity_x_300;

sfwh_velocity_y_10 = open('sfwh_velocity_y_10.dat');
sfwh_velocity_y_10 = sfwh_velocity_y_10.sfwh_velocity_y_10;
sfwh_velocity_y_2 = open('sfwh_velocity_y_2.dat');
sfwh_velocity_y_2 = sfwh_velocity_y_2.sfwh_velocity_y_2;
sfwh_velocity_y_50 = open('sfwh_velocity_y_50.dat');
sfwh_velocity_y_50 = sfwh_velocity_y_50.sfwh_velocity_y_50;
sfwh_velocity_y_100 = open('sfwh_velocity_y_100.dat');
sfwh_velocity_y_100 = sfwh_velocity_y_100.sfwh_velocity_y_100;
sfwh_velocity_y_300 = open('sfwh_velocity_y_300.dat');
sfwh_velocity_y_300 = sfwh_velocity_y_300.sfwh_velocity_y_300;

sfwh_matrix_10(:,:,1) = sfwh_pressure_10;
sfwh_matrix_10(:,:,2) = sfwh_velocity_x_10;
sfwh_matrix_10(:,:,3) = sfwh_velocity_y_10;

sfwh_matrix_2(:,:,1) = sfwh_pressure_2;
sfwh_matrix_2(:,:,2) = sfwh_velocity_x_2;
sfwh_matrix_2(:,:,3) = sfwh_velocity_y_2;

sfwh_matrix_50(:,:,1) = sfwh_pressure_50;
sfwh_matrix_50(:,:,2) = sfwh_velocity_x_50;
sfwh_matrix_50(:,:,3) = sfwh_velocity_y_50;

sfwh_matrix_100(:,:,1) = sfwh_pressure_100;
sfwh_matrix_100(:,:,2) = sfwh_velocity_x_100;
sfwh_matrix_100(:,:,3) = sfwh_velocity_y_100;

sfwh_matrix_300(:,:,1) = sfwh_pressure_300;
sfwh_matrix_300(:,:,2) = sfwh_velocity_x_300;
sfwh_matrix_300(:,:,3) = sfwh_velocity_y_300;

%sfwh_velocity_x = open('sfwh_velocity_x.dat');
%sfwh_velocity_x = sfwh_velocity_x.sfwh_velocity_x;
%sfwh_velocity_y = open('sfwh_velocity_y.dat');
%sfwh_velocity_y = sfwh_velocity_y.sfwh_velocity_y;
% prepare matrix ffowcs-williams H. surface
%sfwh_matrix(:,:,1) = sfwh_pressure;
%sfwh_matrix(:,:,2) = sfwh_velocity_x;
%sfwh_matrix(:,:,3) = sfwh_velocity_y;

% open benchmarks data
DNS_vortex_shedding = open('DNS_vortex_shedding.csv');
DNS_vortex_shedding = DNS_vortex_shedding.DNS_vortex_shedding;
%DNS_vortex_shedding(:,1) = DNS_vortex_shedding(:,1)/max(DNS_vortex_shedding(:,1)); 
%Curle_vortex_shedding = open('Curle_vortex_shedding.csv');
%Curle_vortex_shedding = Curle_vortex_shedding.Curle_vortex_shedding;
%Curle_vortex_shedding(:,1) = Curle_vortex_shedding(:,1)/max(Curle_vortex_shedding(:,1));

f = 0.006; % frequency peak
%theta = (0:10:360)*pi/180;
theta = DNS_vortex_shedding(:,2)*pi/180; % theta in radians
delta_x = 0.004177e-3;
size_square = 2;
radius = 75*(size_square/delta_x);
%radius = 75*2;
pressures_theta_10(1:length(theta)) = 0;
pressures_theta_2(1:length(theta)) = 0;
pressures_theta_50(1:length(theta)) = 0;
pressures_theta_100(1:length(theta)) = 0;
pressures_theta_300(1:length(theta)) = 0;
U = 0.2/sqrt(3);
V = 0;

for point = 1:length(theta)
	disp((point/length(theta)*100));
	[x y] = pol2cart(theta(point), radius);
	pesc(1) = x;
	pesc(2) = y;
	pressures_theta_10(point) = fwh(sfwh_matrix_10,f,U,V,pesc);
	pressures_theta_2(point) = fwh(sfwh_matrix_2,f,U,V,pesc);
	pressures_theta_50(point) = fwh(sfwh_matrix_50,f,U,V,pesc);
	pressures_theta_100(point) = fwh(sfwh_matrix_100,f,U,V,pesc);
	pressures_theta_300(point) = fwh(sfwh_matrix_300,f,U,V,pesc);
end
theta = theta';
pressures_theta_10 = (abs(pressures_theta_10)/(U^2));
%pressures_theta_10 = pressures_theta_10/max(pressures_theta_10);
pressures_theta_10 = pressures_theta_10';
pressures_theta_2 = (abs(pressures_theta_2)/(U^2));
%pressures_theta_2 = pressures_theta_2/max(pressures_theta_10);
pressures_theta_2 = pressures_theta_2';
pressures_theta_50 = (abs(pressures_theta_50)/(U^2));
%pressures_theta_50 = pressures_theta_50/max(pressures_theta_10);
pressures_theta_50 = pressures_theta_50';
pressures_theta_100 = (abs(pressures_theta_100)/(U^2));
%pressures_theta_100 = pressures_theta_100/max(pressures_theta_10);
pressures_theta_100 = pressures_theta_100';
pressures_theta_300 = (abs(pressures_theta_300)/(U^2));
%pressures_theta_300 = pressures_theta_300/max(pressures_theta_10);
pressures_theta_300 = pressures_theta_300';


%compare_pressure_theta_10 = open('pressures_theta.mat');
%compare_pressure_theta_10 = compare_pressure_theta_10.pressures_theta;

sfwh_pressure_5 = open('sfwh_pressure_5.dat');
sfwh_pressure_5 = sfwh_pressure_5.sfwh_pressure_5;
sfwh_velocity_x_5 = open('sfwh_velocity_x_5.dat');
sfwh_velocity_x_5 = sfwh_velocity_x_5.sfwh_velocity_x_5;
sfwh_velocity_y_5 = open('sfwh_velocity_y_5.dat');
sfwh_velocity_y_5 = sfwh_velocity_y_5.sfwh_velocity_y_5;
sfwh_matrix_5(:,:,1) = sfwh_pressure_5;
sfwh_matrix_5(:,:,2) = sfwh_velocity_x_5;
sfwh_matrix_5(:,:,3) = sfwh_velocity_y_5;
pressures_theta_5(1:length(theta)) = 0;
for point = 1:length(theta)
	disp((point/length(theta)*100));
	[x y] = pol2cart(theta(point), radius);
	pesc(1) = x;
	pesc(2) = y;
	pressures_theta_5(point) = fwh(sfwh_matrix_5,f,U,V,pesc);
end
pressures_theta_5 = (abs(pressures_theta_5)/(U^2));
%pressures_theta_5 = pressures_theta_5/max(pressures_theta_10);
pressures_theta_5 = pressures_theta_5';

figure; 
polar(theta, pressures_theta_10/max(pressures_theta_300), 'blue');hold on;
polar(theta, pressures_theta_2/max(pressures_theta_300), 'red'); hold on;
hold on; polar(theta, pressures_theta_50/max(pressures_theta_300), 'black');
hold on; polar(theta, pressures_theta_100/max(pressures_theta_300), 'magenta');
hold on; polar(theta, pressures_theta_300/max(pressures_theta_300), 'cyan');
hold on; polar(theta, pressures_theta_5/max(pressures_theta_300), 'green');
title( ... 
'Directivity with several sizes of FW-H Surface', ... 
'Interpreter','latex','FontSize',16);
xlabel('Angle [Degree]','Interpreter','latex','FontSize',16); 
ylabel('Normalized Magnitude with Surface D = 300','Interpreter','latex','FontSize',16);
legend('D = 10','D = 2', 'D = 50', 'D = 100', 'D = 300', 'D = 5');


%save(filename,variables);

%figure; polar(theta, pressures_theta, '--');
%hold on; polar(theta, DNS_vortex_shedding(:,1), 'red');
%polar(Curle_vortex_shedding(:,2)*(pi/180), Curle_vortex_shedding(:,1), 'black');
%title( ... 
%'Directivity', ... 
%'Interpreter','latex','FontSize',16);
%xlabel('Angle [Degree]','Interpreter','latex','FontSize',16); 
%ylabel('Normalized Magnitude','Interpreter','latex','FontSize',16);
%legend('LBM and FW-H Analogy - Radius Surface = 2','DNS', 'DNS and Curle Analogy');

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