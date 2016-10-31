% source sound case
close all; clear('all'); clc;

addpath(genpath('../'));

% To probes in a boca
history_pressures_file_20 = 'r_20/history_pressures_boca.dat';
history_velocities_file_20 = 'r_20/history_velocities_boca.dat';
radius_20 = 20;
parameters_20 = get_duct_radiation_parameters_boca(history_pressures_file_20, history_velocities_file_20, radius_20);

history_pressures_file_20_fix = 'r_20_fix/history_pressures_boca.dat';
history_velocities_file_20_fix = 'r_20_fix/history_velocities_boca.dat';
radius_20_fix = 20;
parameters_20_fix = get_duct_radiation_parameters_boca(history_pressures_file_20_fix, history_velocities_file_20_fix, radius_20_fix);

% To probes in a 3R
history_pressures_file_20_fix_3r = 'r_20_fix/history_pressures_3r.dat';
history_velocities_file_20_fix_3r = 'r_20_fix/history_velocities_3r.dat';
correction = 0;
position_20_fix_3r = 3*radius_20_fix  + correction;
parameters_20_fix_3r = get_duct_radiation_parameters(history_pressures_file_20_fix_3r, history_velocities_file_20_fix_3r, radius_20_fix, position_20_fix_3r);

history_pressures_file_20_3r = 'r_20/history_pressures_3r.dat';
history_velocities_file_20_3r = 'r_20/history_velocities_3r.dat';
correction = 0;
position_20_3r = 3*radius_20  + correction;
parameters_20_3r = get_duct_radiation_parameters(history_pressures_file_20_3r, history_velocities_file_20_3r, radius_20, position_20_3r);

% To probes in a 6R
history_pressures_file_20_fix_6r = 'r_20_fix/history_pressures_6r.dat';
history_velocities_file_20_fix_6r = 'r_20_fix/history_velocities_6r.dat';
correction_1 = 0;
position_20_fix_6r = 3*radius_20_fix  + correction_1;
parameters_20_fix_6r = get_duct_radiation_parameters(history_pressures_file_20_fix_6r, history_velocities_file_20_fix_6r, radius_20_fix, position_20_fix_6r);

history_pressures_file_20_6r = 'r_20/history_pressures_6r.dat';
history_velocities_file_20_6r = 'r_20/history_velocities_6r.dat';
correction_1 = 0;
position_20_6r = 3*radius_20  + correction_1;
parameters_20_6r = get_duct_radiation_parameters(history_pressures_file_20_6r, history_velocities_file_20_6r, radius_20, position_20_6r);

% -----------------------------------------------------------------
% PLOTING
colours = distinguishable_colors(10);
color_analytical = colours(3:3, 1:3);
color_boca = colours(9:9, 1:3);
color_3r = colours(5:5, 1:3);
color_6r = colours(6:6, 1:3);

% Ploting impedance
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/levine.fig',1);
hold on;
plot(parameters_20_fix{1}, real(parameters_20_fix{2}),'Color' ,[color_boca]); hold on; plot(parameters_20_fix{1}, imag(parameters_20_fix{2}), 'Color' ,[color_boca]);
hold on;
plot(parameters_20_fix_3r{1}, real(parameters_20_fix_3r{2}),'Color' ,[color_3r]); hold on; plot(parameters_20_fix_3r{1}, imag(parameters_20_fix_3r{2}), 'Color' ,[color_3r]);
hold on;
plot(parameters_20_fix_6r{1}, real(parameters_20_fix_6r{2}),'Color' ,[color_6r]); hold on; plot(parameters_20_fix_6r{1}, imag(parameters_20_fix_6r{2}), 'Color' ,[color_6r]);
axis([0 2.5 -0.4 1.2]);
title('Impedance');
ylabel('Impedance','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
legend('Analytical real','Analytical imaginary',...
	' boca real',' boca imaginary', ...
	'3r real','3r imaginary', ...
	'6r real','6r imaginary');
xlim([0 3])
hold off;

% Ploting coefficient of reflection
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/abs_r.fig',1)
hold on;
plot(parameters_20_fix{1}, abs(parameters_20_fix{3}), 'Color', [color_boca]);
hold on;
plot(parameters_20_fix_3r{1}, abs(parameters_20_fix_3r{3}), 'Color', [color_3r]);
hold on;
plot(parameters_20_fix_6r{1}, abs(parameters_20_fix_6r{3}), 'Color', [color_6r]);
%plot(ka,abs(Rr),'--black');
axis([0 2.5 0 1.5]);
title('coefficient of reflection');
ylabel('Coeficiente de Reflexao, Rr','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
legend('Analytical',' boca','3r','6r');
hold off;

% Ploting correction of termination
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/loa.fig',1)
hold on;
plot(parameters_20_fix{1},parameters_20_fix{4}, 'Color', [color_boca]);
hold on;
plot(parameters_20_fix_3r{1},parameters_20_fix_3r{4}, 'Color', [color_3r]);
hold on;
plot(parameters_20_fix_6r{1},parameters_20_fix_6r{4}, 'Color', [color_6r]);
title('correction of termination');
ylabel('End correction, l/a','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
axis([0 1.5 -1 1.2]);
legend('Analytical',' boca','3r','6r');
%legend('Analitico','LBM');

% ----------------------------------------------------------------------
