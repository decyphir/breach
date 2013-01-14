function [mu, mustar, sigm] = EEffects(Y, D, p)
%
% EEFFECTS computes the elementary effects for the global sensitivity method
%          of Morris, given values Y computed at points in a trajectory generated by
%          EE_traj
%
% Synopsis:  [mu mustar sigm] = EEffects(Y, D, p)
%
%  Example:
%
%  P0 = CreateParamSet(Sys,indX);
%  Pr = pRefine(P,4,10);
%  Pf = ComputeTraj(Sys, Pr, tspan);
%  Y = GetTrajValues(Pf, indx, Te);
%
%  [mu mustar sigm] = EEffects(Y, Pr.D, p)


n = size(D,1);
delta = p/(2*(p-1));

% computes diff
dY = diff(Y);

idx = 1:size(Y,2);

% size should be some multiple k of n+1 whereas D is some multiple of n. We
% need to remove indices n+1, 2(n+1), ..., (k-1)*n-1 in dY
idx = idx(mod(idx,n+1)~=0); % abracadabra

dY = dY(idx);

DY = repmat(dY/delta,[n 1]);

EE = DY.*D;
r = numel(Y)/n;

mu = 1/r*sum(EE, 2);
mustar = 1/r* sum(abs(EE), 2);

MU = repmat(mu,[1 size(EE,2)]);
sigm = sqrt(1/(r-1)* sum((EE-MU).^2,2));
