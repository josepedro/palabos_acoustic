% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta função calcula a pressão acústica P em campo distante através da
% superfície bi-dimenssional de Fowcs-Williams e Hawkings, proposta por 
% Parametros a definir Lockard (2000). Os dados de entrada são os
% segu8intes:
%
% P = fwh(MM,f,U,V,pesc)
%
% MM é a matriz obtida através da simulação computacional com os dados de
% pressão total p, velocidade total ux e uy capturados na superfície para
% cada timestep. A matriz MM possui as seguintes dimensões:
% [n_pontos,n_timesteps,3] (pressao, ux e uy). As linhas são os resultados obtidos para cada
% ponto na superfície para um único timestep. As duas primeiras colunas são
% as coordenadas de cada ponto da superfície (x e y). As demais colunas são
% os resultados obtidos na superfície para cada timestep. A primeira folha
% da matriz apresenta a pressão total p. As folhas subsequentes são as
% velocidades totais ux e uy, respectivamente.
%
% f é a frequencia de análise em unidade de lattice (varia de 0 a 0.5)
% U é velocidade média do escoamento em x              
% V é velocidade média do escoamento em y                 
% pesc é ponto de escuta em unidade de lattice (deve ser um vetos do tipo
% [coord_x coord_y]
%
% Andry R. da Silva. UFSC 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = fwh(MM,f,U,V,pesc)


%% Aqui é onde tudo começa

% Propriedades do escoamento


omega = 2*pi*f; % frequ em rad/sec em lattice
k = omega*sqrt(3); % numero de onda


theta = atan(V/U);
M = sqrt(U^2+V^2)*sqrt(3);
beta = sqrt(1-M^2);


x_ = (pesc(1)-MM(:,1,1))*cos(theta) + (pesc(2)-MM(:,2,1))*sin(theta);
y_ = -(pesc(1)-MM(:,1,1))*sin(theta) + (pesc(2)-MM(:,2,1))*cos(theta);

% calculando a funcao de Hankel de segundo tipo e ordem 0
argu = (k/(beta^2))*sqrt(x_.^2+beta^2*y_.^2);
H = besselj(0,argu)-i*bessely(0,argu);

% Calculando a função de Green

G = (i/(4*beta))*exp(i*M*k*x_/(beta^2)).*H;

%% Calculo dos termos monopolo e dipolo

% termo monopolo e F1 e F2
tam = size(MM);
rho_0 = 1;

for col = 3:tam(2)
 
u1 = MM(1:tam(1)-1,col,2);
u2 = MM(1:tam(1)-1,col,3);
p = MM(1:tam(1)-1,col,1);
ro = p*3;
difx = diff(MM(:,1,1));
dify = diff(MM(:,2,1));

% termo monopolo
Q(:,col-2) = (ro.*u1 - rho_0*U).*difx + (ro.*u2 - rho_0*V).*dify;


% termos dipolo
F1(:,col-2) = (p+ro.*(u1-2*U).*u1+rho_0*U^2).*difx + ...
    (ro.*(u1-2*U).*u2 +rho_0*U*V).*dify;

F2(:,col-2) = (ro.*(u2-2*V).*u1 + rho_0*V*U).*difx + ...
    (p+ro.*(u2-2*V).*u2 + rho_0*V^2).*dify;

end

% Fazendo a transformada de fourier de Q e Fs

Q = fft(Q')';
F1 = fft(F1')';
F2 = fft(F2')';


freq = linspace(0,1,tam(2)-1);
[c index] = min(abs(freq-f));
% Calculando as integrais

dGdx = diff(G)./(diff(MM(:,1,1))+eps);
dGdy = diff(G)./(diff(MM(:,2,1))+eps);

P = -sum(F1(:,index).*dGdx) - sum(F2(:,index).*dGdy) - sum(i*omega*Q(:,index).*G(1:length(G)-1));





