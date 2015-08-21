clc;
clearvars;
close all;

tic;                    % calculate program run tim

% helper functions to convert between frequency and wavelength
c = 3e8;
f2w = @(freq) c/freq;
w2f = @(wave) c/wave;
w2f10 = @(wave) log10(c/wave);

dir = 'linear';
N =21;
wave = 0.59e-6;
freq = w2f(wave);    	% operating wavelength
mu_r = 1;               % relative permeability
epi_r = 1;              % relative permittivity

% incident e field
E_inc = 1;             
E_vect = [0 1 0];
% E_vect = [sin((90-120/2)*2*pi/360) cos((90-120/2)*2*pi/360) 0];  
% assump = 0;             % assume end currents are zero
assump = 1;             % assume end currents aren't zero

% analyze antenna setup
ant1 = Antenna(dir, N, freq, epi_r, mu_r, E_inc, E_vect, assump);
% ant1 = Antenna_noends(dir, N, freq, epi_r, mu_r, E_inc, E_vect);

toc

fid = fopen('blah.csv', 'wt');
fprintf(fid, '%s, %e, %e', dir, L,a);
% xz plane 
% xIn = -1.5e-6:0.1e-6:1.5e-6;
% yIn = 0;
% zIn = -1e-6:0.1e-6:9e-6;

% xy plane
% xIn = -2.5e-6:0.1e-6:2.5e-6;
% yIn = -2.5e-6:0.1e-6:2.5e-6;
% zIn = 0e-6;

xIn = -6e-6:0.2e-6:6e-6;
yIn = -6e-6:0.2e-6:6e-6;
zIn = 1e-6;
% 
% ant1.plotField(xIn, yIn, zIn, 1);
