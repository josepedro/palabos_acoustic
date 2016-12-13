% Analysing coefficient of reflection

close all; clear('all'); clc;

addpath(genpath('../'));

length_buffer_anechoic = 30;

% To probes
pressures_point_A = textread('point_A/history_pressures_point_A.dat');
velocities_x_point_A = 0;
velocities_y_point_A = 0;
velocities_z_point_A = textread('point_A/history_velocities_z_point_A.dat');
parameters_point_A = get_anechoic_parameters(pressures_point_A, velocities_x_point_A, velocities_y_point_A, velocities_z_point_A, length_buffer_anechoic);

pressures_point_B = textread('point_B/history_pressures_point_B.dat');
velocities_x_point_B = 0;
velocities_y_point_B = textread('point_B/history_velocities_y_point_B.dat');
velocities_z_point_B = textread('point_B/history_velocities_z_point_B.dat');
parameters_point_B = get_anechoic_parameters(pressures_point_B, velocities_x_point_B, velocities_y_point_B, velocities_z_point_B, length_buffer_anechoic);

pressures_point_C = textread('point_C/history_pressures_point_C.dat');
velocities_x_point_C = 0;
velocities_y_point_C = textread('point_C/history_velocities_y_point_C.dat');
velocities_z_point_C = 0;
parameters_point_C = get_anechoic_parameters(pressures_point_C, velocities_x_point_C, velocities_y_point_C, velocities_y_point_C, length_buffer_anechoic);

pressures_point_D = textread('point_D/history_pressures_point_D.dat');
velocities_x_point_D = 0;
velocities_y_point_D = 0;
velocities_z_point_D = -textread('point_D/history_velocities_z_point_D.dat');
parameters_point_D = get_anechoic_parameters(pressures_point_D, velocities_x_point_D, velocities_y_point_D, velocities_z_point_D, length_buffer_anechoic);

% -----------------------------------------------------------------
% PLOTING
colours = distinguishable_colors(4);
color_point_A = colours(1:1, 1:3);
color_point_B = colours(2:2, 1:3);
color_point_C = colours(3:3, 1:3);
color_point_D = colours(4:4, 1:3);

% Ploting impedance
figure(1);
hold on;
plot(parameters_point_A{1}, real(parameters_point_A{2}),'Color' ,[color_point_A]); hold on; plot(parameters_point_A{1}, imag(parameters_point_A{2}), '--', 'Color' ,[color_point_A]);
hold on;
plot(parameters_point_B{1}, real(parameters_point_B{2}),'Color' ,[color_point_B]); hold on; plot(parameters_point_B{1}, imag(parameters_point_B{2}), '--', 'Color' ,[color_point_B]);
hold on;
plot(parameters_point_C{1}, real(parameters_point_C{2}),'Color' ,[color_point_C]); hold on; plot(parameters_point_C{1}, imag(parameters_point_C{2}), '--', 'Color' ,[color_point_C]);
hold on;
plot(parameters_point_D{1}, real(parameters_point_D{2}),'Color' ,[color_point_D]); hold on; plot(parameters_point_D{1}, imag(parameters_point_D{2}), '--', 'Color' ,[color_point_D]);
axis([0 2.5 -0.4 1.2]);
title('Impedance');
ylabel('Impedance','FontSize',20);
xlabel('Numero de Onda, kl','FontSize',20);
legend('ponto A real','ponto A imaginario',...
	'ponto B real','ponto B imaginario',...
	'ponto C real','ponto C imaginario',...
	'ponto D real','ponto D imaginario');
xlim([0 3])
hold off;

% Ploting coefficient of reflection
figure(2);
hold on;
plot(parameters_point_A{1}, abs(parameters_point_A{3}), 'Color', [color_point_A]);
hold on;
plot(parameters_point_B{1}, abs(parameters_point_B{3}), 'Color', [color_point_B]);
hold on;
plot(parameters_point_C{1}, abs(parameters_point_C{3}), 'Color', [color_point_C]);
hold on;
plot(parameters_point_D{1}, abs(parameters_point_D{3}), 'Color', [color_point_D]);
%plot(ka,abs(Rr),'--black');
%axis([0 2.5 0 1.5]);
xlim([0 2.5])
title('coefficient of reflection');
ylabel('Coeficiente de Reflexao, Rr','FontSize',20);
xlabel('Numero de Onda, ka','FontSize',20);
legend('ponto A','ponto B','ponto C','ponto D');
hold off;


% ----------------------------------------------------------------------
% SAVING FIGURES

saveas(1, 'impedance', 'png');
saveas(1, 'impedance', 'fig');
saveas(2, 'abs_r', 'png');
saveas(2, 'abs_r', 'fig');