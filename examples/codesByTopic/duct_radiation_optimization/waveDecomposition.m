function [pplus,pminus] = waveDecomposition(pressure,positions,refPos,kplus,kminus)
	% PROPAGATINGWAVES is the implementation of equation 4 to find the waves
	% propagating in both directions using a pseudo-inverse approach from
	% pressure measurements at multiple locations.
	% r = radius;
	% press=101325;
	% c0=331*sqrt(temperature+273)/16.52;            % Speed of Sound
	% rho0=press/(287*(temperature+273));            % Density of air
	% Pr=0.7;                                 %Prandtl Number
	% Tk=temperature+273;
	% V1 = 0.784812; V2 = -6.75085E-04; v3 = 4.11244E-07; v4 = -1.0343E-10;
	% mu=(Tk.*(V1+Tk.*(V2+Tk.*(v3+Tk*v4))))*0.0000001;
	% gamma = 1.4022; %specific heats ratio
	% 
	% w = frequency*2*pi;
	% s = r*sqrt(rho0*w/mu); %shear wavenumber
	% K0 = 1+ ((1-1j)/(s*sqrt(2)))*(1+(gamma-1)/sqrt(Pr));
	% k0 = w/c0;
	% kplus = k0*K0/(1+K0*M);
	% kminus = k0*(K0)/(1-K0*M);

	e = zeros(length(positions),2);
	for n = 1:length(positions);
	    e(n,1) = exp(-1i*kplus*(positions(n)-refPos)); % P(1) = a(1,1) * P(2)
	    e(n,2) = exp(1i*kminus*(positions(n)-refPos)); %no wave propagation in z-
	end
	% Least Squares approach using MP pseudo-inverse
	Pfinal = (pinv(e)*pressure);
	pplus = Pfinal(1);
	pminus = Pfinal(2);
end