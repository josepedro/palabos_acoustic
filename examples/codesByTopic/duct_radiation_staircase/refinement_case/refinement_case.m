% refinement case
close all; clear('all'); clc;

addpath(genpath('../'));

% To probes in a boca
history_pressures_file_10 = 'r_10/history_pressures_boca.dat';
history_velocities_file_10 = 'r_10/history_velocities_boca.dat';
radius_10 = 10;
parameters_10 = get_duct_radiation_parameters_boca(history_pressures_file_10, history_velocities_file_10, radius_10);

history_pressures_file_20 = 'r_20/history_pressures_boca.dat';
history_velocities_file_20 = 'r_20/history_velocities_boca.dat';
radius_20 = 20;
parameters_20 = get_duct_radiation_parameters_boca(history_pressures_file_20, history_velocities_file_20, radius_20);

history_pressures_file_30 = 'r_30/history_pressures_boca.dat';
history_velocities_file_30 = 'r_30/history_velocities_boca.dat';
radius_30 = 30;
parameters_30 = get_duct_radiation_parameters_boca(history_pressures_file_30, history_velocities_file_30, radius_30);

history_pressures_file_40 = 'r_40/history_pressures_boca.dat';
history_velocities_file_40 = 'r_40/history_velocities_boca.dat';
radius_40 = 40;
parameters_40 = get_duct_radiation_parameters_boca(history_pressures_file_40, history_velocities_file_40, radius_40);

% Ploting impedance
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/levine.fig',1);
hold on;
plot(parameters_10{1}, real(parameters_10{2}),'--'); hold on; plot(parameters_10{1}, imag(parameters_10{2}), 'r--');
hold on;
plot(parameters_20{1}, real(parameters_20{2}),'o'); hold on; plot(parameters_20{1}, imag(parameters_20{2}), 'ro');
hold on;
plot(parameters_30{1}, real(parameters_30{2}),'*'); hold on; plot(parameters_30{1}, imag(parameters_30{2}), 'r*');
hold on;
plot(parameters_40{1}, real(parameters_40{2}),'+'); hold on; plot(parameters_40{1}, imag(parameters_40{2}), 'r+');
%axis([0 1.82 0 max(abs(ZL))]);
title('Impedance boca');
ylabel('Impedance','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
%legend('Impedance real part','Impedance imaginary part');
xlim([0 3])
hold off;

% Ploting coefficient of reflection
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/abs_r.fig',1)
hold on;
plot(parameters_10{1}, abs(parameters_10{3}), '--');
hold on;
plot(parameters_20{1}, abs(parameters_20{3}), 'o');
hold on;
plot(parameters_30{1}, abs(parameters_30{3}), '*');
hold on;
plot(parameters_40{1}, abs(parameters_40{3}), '+');
%plot(ka,abs(Rr),'--black');
axis([0 2.5 0 1.5]);
title('coefficient of reflection boca');
ylabel('Coeficiente de Reflexao, Rr','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
%legend('Analitico','LBM');
hold off;

% Ploting correction of termination
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/loa.fig',1)
hold on;
plot(parameters_10{1},parameters_10{4}, '--');
hold on;
plot(parameters_20{1},parameters_20{4}, 'o');
hold on;
plot(parameters_30{1},parameters_30{4}, '*');
hold on;
plot(parameters_40{1},parameters_40{4}, '+');
title('coefficient of reflection boca');
ylabel('End correction, l/a','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
axis([0 1.5 -1 1]);
%legend('Analitico','LBM');


% To probes in a 3R
history_pressures_file_10_3r = 'r_10/history_pressures_3r.dat';
history_velocities_file_10_3r = 'r_10/history_velocities_3r.dat';
correction = 0;
position_10_3r = 3*radius_10  + correction;
parameters_10_3r = get_duct_radiation_parameters(history_pressures_file_10_3r, history_velocities_file_10_3r, radius_10, position_10_3r);

history_pressures_file_20_3r = 'r_20/history_pressures_3r.dat';
history_velocities_file_20_3r = 'r_20/history_velocities_3r.dat';
position_20_3r = 3*radius_20  + correction;
parameters_20_3r = get_duct_radiation_parameters(history_pressures_file_20_3r, history_velocities_file_20_3r, radius_20, position_20_3r);

history_pressures_file_30_3r = 'r_30/history_pressures_3r.dat';
history_velocities_file_30_3r = 'r_30/history_velocities_3r.dat';
position_30_3r = 3*radius_30  + correction;
parameters_30_3r = get_duct_radiation_parameters(history_pressures_file_30_3r, history_velocities_file_30_3r, radius_30, position_30_3r);

history_pressures_file_40_3r = 'r_40/history_pressures_3r.dat';
history_velocities_file_40_3r = 'r_40/history_velocities_3r.dat';
position_40_3r = 3*radius_40  + correction;
parameters_40_3r = get_duct_radiation_parameters(history_pressures_file_40_3r, history_velocities_file_40_3r, radius_40, position_40_3r);

% Ploting impedance
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/levine.fig',1);
hold on;
plot(parameters_10_3r{1}, real(parameters_10_3r{2}),'--'); hold on; plot(parameters_10_3r{1}, imag(parameters_10_3r{2}), 'r--');
hold on;
plot(parameters_20_3r{1}, real(parameters_20_3r{2}),'o'); hold on; plot(parameters_20_3r{1}, imag(parameters_20_3r{2}), 'ro');
hold on;
plot(parameters_30_3r{1}, real(parameters_30_3r{2}),'*'); hold on; plot(parameters_30_3r{1}, imag(parameters_30_3r{2}), 'r*');
hold on;
plot(parameters_40_3r{1}, real(parameters_40_3r{2}),'+'); hold on; plot(parameters_40_3r{1}, imag(parameters_40_3r{2}), 'r+');
%axis([0 1.82 0 max(abs(ZL))]);
ylabel('Impedance','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
%legend('Impedance real part','Impedance imaginary part');
title('Impedance 3r');
xlim([0 3])
hold off;

% Ploting coefficient of reflection
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/abs_r.fig',1)
hold on;
plot(parameters_10_3r{1}, abs(parameters_10_3r{3}), '--');
hold on;
plot(parameters_20_3r{1}, abs(parameters_20_3r{3}), 'o');
hold on;
plot(parameters_30_3r{1}, abs(parameters_30_3r{3}), '*');
hold on;
plot(parameters_40_3r{1}, abs(parameters_40_3r{3}), '+');
%plot(ka,abs(Rr),'--black');
axis([0 2.5 0 1.5]);
title('coefficient of reflection 3r');
ylabel('Coeficiente de Reflexao, Rr','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
%legend('Analitico','LBM');
hold off;

% Ploting correction of termination
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/loa.fig',1)
hold on;
plot(parameters_10_3r{1},parameters_10_3r{4}, '--');
hold on;
plot(parameters_20_3r{1},parameters_20_3r{4}, 'o');
hold on;
plot(parameters_30_3r{1},parameters_30_3r{4}, '*');
hold on;
plot(parameters_40_3r{1},parameters_40_3r{4}, '+');
title('correction of termination 3r');
ylabel('End correction, l/a','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
axis([0 1.5 -1 1]);
%legend('Analitico','LBM');


% To probes in a 6R
history_pressures_file_10_6r = 'r_10/history_pressures_6r.dat';
history_velocities_file_10_6r = 'r_10/history_velocities_6r.dat';
correction_1 = 0;
position_10_6r = 6*radius_10  + correction_1;
parameters_10_6r = get_duct_radiation_parameters(history_pressures_file_10_6r, history_velocities_file_10_6r, radius_10, position_10_6r);

history_pressures_file_20_6r = 'r_20/history_pressures_6r.dat';
history_velocities_file_20_6r = 'r_20/history_velocities_6r.dat';
position_20_6r = 3*radius_20  + correction_1;
parameters_20_6r = get_duct_radiation_parameters(history_pressures_file_20_6r, history_velocities_file_20_6r, radius_20, position_20_6r);

history_pressures_file_30_6r = 'r_30/history_pressures_6r.dat';
history_velocities_file_30_6r = 'r_30/history_velocities_6r.dat';
position_30_6r = 3*radius_30  + correction_1;
parameters_30_6r = get_duct_radiation_parameters(history_pressures_file_30_6r, history_velocities_file_30_6r, radius_30, position_30_6r);

history_pressures_file_40_6r = 'r_40/history_pressures_6r.dat';
history_velocities_file_40_6r = 'r_40/history_velocities_6r.dat';
position_40_6r = 3*radius_40  + correction_1;
parameters_40_6r = get_duct_radiation_parameters(history_pressures_file_40_6r, history_velocities_file_40_6r, radius_40, position_40_6r);

% Ploting impedance
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/levine.fig',1);
hold on;
plot(parameters_10_6r{1}, real(parameters_10_6r{2}),'--'); hold on; plot(parameters_10_6r{1}, imag(parameters_10_6r{2}), 'r--');
hold on;
plot(parameters_20_6r{1}, real(parameters_20_6r{2}),'o'); hold on; plot(parameters_20_6r{1}, imag(parameters_20_6r{2}), 'ro');
hold on;
plot(parameters_30_6r{1}, real(parameters_30_6r{2}),'*'); hold on; plot(parameters_30_6r{1}, imag(parameters_30_6r{2}), 'r*');
hold on;
plot(parameters_40_6r{1}, real(parameters_40_6r{2}),'+'); hold on; plot(parameters_40_6r{1}, imag(parameters_40_6r{2}), 'r+');
%axis([0 1.82 0 max(abs(ZL))]);
ylabel('Impedance','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
%legend('Impedance real part','Impedance imaginary part');
title('Impedance 6r');
xlim([0 3])
hold off;

% Ploting coefficient of reflection
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/abs_r.fig',1)
hold on;
plot(parameters_10_6r{1}, abs(parameters_10_6r{3}), '--');
hold on;
plot(parameters_20_6r{1}, abs(parameters_20_6r{3}), 'o');
hold on;
plot(parameters_30_6r{1}, abs(parameters_30_6r{3}), '*');
hold on;
plot(parameters_40_6r{1}, abs(parameters_40_6r{3}), '+');
%plot(ka,abs(Rr),'--black');
axis([0 2.5 0 1.5]);
title('coefficient of reflection 6r');
ylabel('Coeficiente de Reflexao, Rr','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
%legend('Analitico','LBM');
hold off;

% Ploting correction of termination
uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/loa.fig',1)
hold on;
plot(parameters_10_6r{1},parameters_10_6r{4}, '--');
hold on;
plot(parameters_20_6r{1},parameters_20_6r{4}, 'o');
hold on;
plot(parameters_30_6r{1},parameters_30_6r{4}, '*');
hold on;
plot(parameters_40_6r{1},parameters_40_6r{4}, '+');
title('correction of termination 6r');
ylabel('End correction, l/a','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
axis([0 1.5 -1 1]);
%legend('Analitico','LBM');

