% Analysing coefficient of reflection

close all; clear('all'); clc;

addpath(genpath('../'));

length_buffer_anechoic = 30;

% To probes
pressures_abc_A = textread('abc_A/history_pressures_abc_A.dat');
velocities_x_abc_A = 0;
velocities_y_abc_A = 0;
velocities_z_abc_A = textread('abc_A/history_velocities_z_abc_A.dat');
parameters_abc_A = get_anechoic_parameters(pressures_abc_A, velocities_x_abc_A, velocities_y_abc_A, velocities_z_abc_A, length_buffer_anechoic);

pressures_abc_B = textread('abc_B/history_pressures_abc_B.dat');
velocities_x_abc_B = 0;
velocities_y_abc_B = textread('abc_B/history_velocities_y_abc_B.dat');
velocities_z_abc_B = textread('abc_B/history_velocities_z_abc_B.dat');
parameters_abc_B = get_anechoic_parameters(pressures_abc_B, velocities_x_abc_B, velocities_y_abc_B, velocities_z_abc_B, length_buffer_anechoic);

pressures_abc_C = textread('abc_C/history_pressures_abc_C.dat');
velocities_x_abc_C = 0;
velocities_y_abc_C = textread('abc_C/history_velocities_y_abc_C.dat');
velocities_z_abc_C = 0;
parameters_abc_C = get_anechoic_parameters(pressures_abc_C, velocities_x_abc_C, velocities_y_abc_C, velocities_y_abc_C, length_buffer_anechoic);

pressures_abc_D = textread('abc_D/history_pressures_abc_D.dat');
velocities_x_abc_D = 0;
velocities_y_abc_D = 0;
velocities_z_abc_D = -textread('abc_D/history_velocities_z_abc_D.dat');
parameters_abc_D = get_anechoic_parameters(pressures_abc_D, velocities_x_abc_D, velocities_y_abc_D, velocities_z_abc_D, length_buffer_anechoic);

% -----------------------------------------------------------------
% PLOTING
colours = distinguishable_colors(4);
color_abc_A = colours(1:1, 1:3);
color_abc_B = colours(2:2, 1:3);
color_abc_C = colours(3:3, 1:3);
color_abc_D = colours(4:4, 1:3);

% Ploting impedance
figure(1);
hold on;
plot(parameters_abc_A{1}, real(parameters_abc_A{2}),'Color' ,[color_abc_A]); hold on; plot(parameters_abc_A{1}, imag(parameters_abc_A{2}), '--', 'Color' ,[color_abc_A]);
hold on;
plot(parameters_abc_B{1}, real(parameters_abc_B{2}),'Color' ,[color_abc_B]); hold on; plot(parameters_abc_B{1}, imag(parameters_abc_B{2}), '--', 'Color' ,[color_abc_B]);
hold on;
plot(parameters_abc_C{1}, real(parameters_abc_C{2}),'Color' ,[color_abc_C]); hold on; plot(parameters_abc_C{1}, imag(parameters_abc_C{2}), '--', 'Color' ,[color_abc_C]);
hold on;
%plot(parameters_abc_D{1}, real(parameters_abc_D{2}),'Color' ,[color_abc_D]); hold on; plot(parameters_abc_D{1}, imag(parameters_abc_D{2}), '--', 'Color' ,[color_abc_D]);
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
plot(parameters_abc_A{1}, abs(parameters_abc_A{3}), 'Color', [color_abc_A]);
hold on;
plot(parameters_abc_B{1}, abs(parameters_abc_B{3}), 'Color', [color_abc_B]);
hold on;
%plot(parameters_abc_C{1}, abs(parameters_abc_C{3}), 'Color', [color_abc_C]);
hold on;
%plot(parameters_abc_D{1}, abs(parameters_abc_D{3}), 'Color', [color_abc_D]);
%plot(ka,abs(Rr),'--black');
%axis([0 2.5 0 1.5]);
xlim([0 0.5])
%ylim([0 2])
title('coefficient of reflection');
ylabel('Coeficiente de Reflexao, Rr','FontSize',20);
xlabel('Frequencias, ka','FontSize',20);
legend('ponto A','ponto B','ponto C','ponto D');
hold off;


% ----------------------------------------------------------------------
% SAVING FIGURES

saveas(1, 'impedance_abc', 'png');
saveas(1, 'impedance_abc', 'fig');
saveas(2, 'abs_r_abc', 'png');
saveas(2, 'abs_r_abc', 'fig');