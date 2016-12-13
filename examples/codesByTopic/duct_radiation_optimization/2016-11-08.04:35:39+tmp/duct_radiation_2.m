% duct_radiation_2

clear('all'); clc; close all;

pressures_1 = textread('2r_p1/history_pressures_2r_p1.dat');
pressures_1 = [pressures_1(:)];
pressures_2 = textread('2r_p2/history_pressures_2r_p2.dat');
pressures_2 = [pressures_2(:)];
type_mic = 2;

window = hanning(length(pressures_1));
pressures_1 = window.*pressures_1;
pressures_2 = window.*pressures_2;

a = 20;
const = 4;
L_1 = type_mic*a + const;
L_2 = type_mic*a - 5  + const;
cs = 1/sqrt(3);
frequencies = linspace(0, 1, length(pressures_1));
frequencies = frequencies';
ka = (2*pi*frequencies*a)/cs;
k = ka/a;
rho0 = 1;

H_12 = (fft(pressures_2))./(fft(pressures_1));
s = L_1 - L_2;
H_I = exp(-i*k*s);
H_R = exp(i*k*s);
A = H_12 - H_I;
B = H_R - H_12;
r = ((A./B)).*exp(2*i*k*L_1);
l = (log(-r./(abs(r))))./(-2*i*k);
l_line = real(l);
la = l_line/(a);
Zr = (1 + r)./(1 - r);


uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/levine.fig',1)
hold on
plot(ka, real(Zr),'--'); hold on; plot(ka, imag(Zr), 'r--');
xlim([0 1.8])
hold off
ylabel('Impedance','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
legend('Impedance real part','Impedance imaginary part');

uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/abs_r.fig',2)
hold on;
%plot(ka,abs(RL),'--black');
plot(ka,abs(r),'--black');
axis([0 2.5 0 1.5]);
ylabel('Coeficiente de Reflex\E3o, Rr','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
%legend('Analitico','LBM');
hold off

uiopen('/home/josepedro/palabos_acoustic/examples/codesByTopic/duct_radiation_staircase/loa.fig',3)
% %load munt_loa.mat
% %plot(loa_005(:,1),loa_005(:,2),'black')
hold on
plot(ka,la)
ylabel('End correction, l/a','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
axis([0 1.5 -1 1]);
legend('Analitico','LBM');