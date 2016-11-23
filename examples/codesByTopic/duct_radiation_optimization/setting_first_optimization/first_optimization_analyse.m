% first optimization analyse

close all; clear('all'); clc;

addpath(genpath('../'));

% To probes in a boca
history_pressures_file_with = 'with/history_pressures_boca.dat';
history_velocities_file_with = 'with/history_velocities_boca.dat';
radius_with = 10;
parameters_with = get_duct_radiation_parameters_boca(history_pressures_file_with, history_velocities_file_with, radius_with);

history_pressures_file_without = 'without/history_pressures_boca.dat';
history_velocities_file_without = 'without/history_velocities_boca.dat';
radius_without = 10;
parameters_without = get_duct_radiation_parameters_boca(history_pressures_file_without, history_velocities_file_without, radius_without);

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
plot(parameters_with{1}, real(parameters_with{2}),'Color' ,[color_boca]); hold on; plot(parameters_with{1}, imag(parameters_with{2}), 'Color' ,[color_boca]);
hold on;
plot(parameters_without{1}, real(parameters_without{2}),'Color' ,[color_3r]); hold on; plot(parameters_without{1}, imag(parameters_without{2}), 'Color' ,[color_3r]);
axis([0 2.5 -0.4 1.2]);
title('Impedance');
ylabel('Impedance','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
legend('Analytical real','Analytical imaginary',...
	' novo real','novo imaginary', ...
	'antigo real','antigo imaginary');
xlim([0 3])
hold off;

% Ploting coefficient of reflection
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/abs_r.fig',1)
hold on;
plot(parameters_with{1}, abs(parameters_with{3}), 'Color', [color_boca]);
hold on;
plot(parameters_without{1}, abs(parameters_without{3}), 'Color', [color_3r]);
%plot(ka,abs(Rr),'--black');
axis([0 2.5 0 1.5]);
title('coefficient of reflection');
ylabel('Coeficiente de Reflexao, Rr','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
legend('Analytical','novo','antigo');
hold off;

% Ploting correction of termination
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/loa.fig',1)
hold on;
plot(parameters_with{1},parameters_with{4}, 'Color', [color_boca]);
hold on;
plot(parameters_without{1},parameters_without{4}, 'Color', [color_3r]);
title('correction of termination');
ylabel('End correction, l/a','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
axis([0 1.5 -1 1.2]);
legend('Analytical','novo','antigo');
%legend('Analitico','LBM');

% ----------------------------------------------------------------------
% SAVING FIGURES

saveas(1, 'impedance', 'png');
saveas(1, 'impedance', 'fig');
saveas(2, 'abs_r', 'png');
saveas(2, 'abs_r', 'fig');
saveas(3, 'loa', 'png');
saveas(3, 'loa', 'fig');
