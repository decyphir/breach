% ************************************************************************
% 2nd Order Sigma-Delta BandPass A/D Modulator
% by S. Brigati & F. Francesconi Ver.(0.1) 24/09/99
% The modulator structure is simulated using Simulink (bp2.mdl).
% Post-processing of the results is done with Matlab.
% 1. Plots the Power Spectral Density of the bit-stream
% 2. Calculates the SNR
% 3. Calculates histograms at the integrator outputs
% ************************************************************************

t0=clock;

% ************************************************************************
% Variabili globali
% ************************************************************************
bw=200e3;				% Base-band
R=107;					% 42.8 MHz sampling rate
% R=50;					% 20.0 MHz sampling rate
Fs=R*2*bw;				% Oversampling frequency
Ts=1/Fs;
N=65536;				% Samples number
nper=4;
Fin=(Fs/4)+ nper*Fs/N;			% Input signal frequency (Fin=(Fs/4)+ nper*Fs/N)
Ampl=0.5-pi/256;		% Input signal amplitude [V]
Ntransient=0;
k=1.38e-23;				% Boltzmann Constant
Temp=300;				% Absolute Temperature in Kelvin
% Temp=0;				% Absolute Temperature in Kelvin
Cf=4e-12;				% Integrating Capacitance of the first integrator
alfa=(711-1)/711;		% A=Op-amp finite gain (alfa=(A-1)/A -> ideal op-amp alfa=1)
% alfa=1;
Amax=2;					% Op-amp saturation value [V]
% Amax=10000;					% Op-amp saturation value [V]
sr=280e6;				% Op-amp slew rate [V/s]
% sr=280e10;				% Op-amp slew rate [V/s]
GBW=250e6;				% Op-amp GBW [Hz]
% GBW=250e10;				% Op-amp GBW [Hz]
noise1=4.39e-3;			% 1st int. output noise std. dev. [Vrms]
% noise1=0;
delta=0.1e-9;				% Random Sampling jitter (std. dev.) [s] (Boser, Wooley JSSC Dec. 88)
% delta=0;
%echo off;
% Modulator coefficients

b=0.125;		% 1/8 gain of the first stage
b2=0.125;	% 1/8 gain of the second stage
b3=0.25;		% 1/4 additional gain of the second feedback loop
Vref=1;
finrad=Fin*2*pi;		% Input signal frequency in radians

% s0=sprintf('** Simulation Parameters **');
% s1=sprintf('   Fs(Hz)=%1.0f',Fs);
% s2=sprintf('   Ts(s)=%1.6e',Ts);
% s3=sprintf('   Fin(Hz)=%1.4f',Fin);
% s4=sprintf('   BW(Hz)=%1.0f',bw);
% s5=sprintf('   OSR=%1.0f',R);
% s6=sprintf('   Npoints=%1.0f',N);
% s7=sprintf('   tsim(sec)=%1.3f',N/Fs);
% s8=sprintf('   Nperiods=%1.3f',N*Fin/Fs);
% disp(s0)
% disp(s1)
% disp(s2)
% disp(s3)
% disp(s4)
% disp(s5)
% disp(s6)
% disp(s7)
% disp(s8)
% ************************************************************************
% Open Simulink diagram first
% ************************************************************************

% options=simset('InitialState', zeros(1,4), 'RelTol', 1e-3, 'MaxStep', 1/Fs);
% sim('bp2', (N+Ntransient)/Fs, options);	% Starts Simulink simulation

% ************************************************************************
%   Calculates SNR and PSD of the bit-stream and of the signal
% ************************************************************************
w=hann(N);
f=Fin/Fs;			% Normalized signal frequency
fB=N*(bw/Fs);		% Base-band frequency bins
fBL=N*(1/4-bw/(2*Fs));		% Lower limit Base-band frequency bins
fBH=N*(1/4+bw/(2*Fs));		% Higher limit Base-band frequency bins