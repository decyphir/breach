%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hsdata.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [u,v,F1,F2,x,fglob,xglob] = hsdata(prob)
% gives the bounds on the variables and constraints, the standard
% starting point, the optimal function value and the global minimizer 
% for the Hock-Schittkowski function number nprob
%
% Input:
% nprob  number of problem
%        at present, nprob in {18,19,23,30,31,34,36,37,41,53,60,65,66,
%        71,73,74} is implemented
% 
% Output:
% u,v    bounds on the variables
% F1,F2  bounds on the constraints
% x      standard starting point
% fglob  optimal objective function value
% xglob  global minimizer
%
function [u,v,F1,F2,x,fglob,xglob] = hsbounds(prob)
if prob == 18
  u = [2; 0];
  v = [50; 50];
  F1 = [0; 0];
  F2 = [Inf; Inf];
  x = [2; 2];
  fglob = 5;
  xglob = [sqrt(250); sqrt(2.5)];
elseif prob == 19
  u = [13; 0];
  v = [100; 100];
  F1 = [0; 0];
  F2 = [Inf; Inf];
  x = [20.1; 5.84];
  fglob = -6961.81381;
  xglob = [14.095; 0.84296079];
elseif prob == 23
  u = [-50; -50];
  v = [50; 50];
  F1 = zeros(5,1);
  F2 = Inf*ones(5,1);
  x = [3; 1];
  fglob  = 2;
  xglob = [1; 1];
elseif prob == 30
  u = [1; -10; -10];
  v = 10*ones(3,1);
  F1 = 0;
  F2 = Inf;
  x = [1; 1; 1];
  fglob = 1;
  xglob = [1; 0; 0];
elseif prob == 31
  u = [-10; 1; -10];
  v = [10; 10; 1];
  F1 = 0;
  F2 = Inf;
  x = [1; 1; 1];
  fglob = 6;
  xglob = [1/sqrt(3); sqrt(3); 0];
elseif prob == 34
  u = zeros(3,1);
  v = [100; 100; 10];
  F1 = [0; 0];
  F2 = [Inf; Inf];
  x = [0; 1.05; 2.9];
  fglob = -log(log(10));
  xglob = [log(log(10)); log(10); 10];
elseif prob == 36
  u = zeros(3,1);
  v = [20; 11; 42];
  F1 = 0;
  F2 = Inf;
  x = [10; 10; 10];
  fglob = -3300;
  xglob = [20; 11; 15];
elseif prob == 37
  u = zeros(3,1);
  v = 42*ones(3,1);
  F1 = [0; 0];
  F2 = [Inf; Inf];
  x = [10; 10; 10];
  fglob = -3456;
  xglob = [24; 12; 12];
elseif prob == 41
  u = zeros(4,1);
  v = [ones(3,1); 2];
  F1 = 0;
  F2 = 0;
  x = [2; 2; 2; 2];
  fglob = 52/27;
  xglob = [2/3; 1/3; 1/3; 2];
elseif prob == 53
  u = -10*ones(5,1);
  v = 10*ones(5,1);
  F1 = zeros(3,1);
  F2 = zeros(3,1);
  x = [2; 2; 2; 2; 2];
  fglob = 176/43;
  xglob = [-33; 11; 27; -5; 11]/43;
elseif prob == 60
  u = -10*ones(3,1);
  v = 10*ones(3,1);
  F1 = 0;
  F2 = 0;
  x = [2; 2; 2];
  fglob = 0.03256820025;
  xglob = [1.104859024; 1.196674194; 1.535262257];
elseif prob == 65
  u = [-4.5; -4.5; -5];
  v = [4.5; 4.5; 5];
  F1 = 0;
  F2 = Inf;
  x = [-5; 5; 0];
  fglob = 0.9535288567;
  xglob = [3.650461821; 3.65046169; 4.6204170507];
elseif prob == 66
  u = zeros(3,1);
  v = [100; 100; 10];
  F1 = [0; 0];
  F2 = [Inf; Inf];
  x = [0; 1.05; 2.9];
  fglob = 0.5181632741;
  xglob = [0.1841264879; 1.202167873; 3.327322322];
elseif prob == 71
  u = ones(4,1);
  v = 5*ones(4,1);
  F1 = [0; 0];
  F2 = [Inf; 0];
  x = [1; 5; 5; 1];
  fglob = 17.0140173;
  xglob = [1; 4.7429994; 3.8211503; 1.3794082];
elseif prob == 73
  u = zeros(4,1);
  v = Inf*ones(4,1);
  F1 = zeros(3,1);
  F2 = [Inf; Inf; 0];
  x = [1; 1; 1; 1];
  fglob = 130.8;
  xglob = [0.6355216; -0.12e-11; 0.3127019; 0.05177655];
elseif prob == 74
  a = 0.55;
  u = [0; 0; -a; -a];
  v = [1200; 1200; a; a];
  F1 = zeros(5,1);
  F2 = [Inf; Inf; zeros(3,1)];
  x = zeros(4,1);
  fglob = 5126.4981;
  xglob = [679.9453; 1026.067; 0.1188764; -0.3962336];
else
  error('Problem not yet coded!')
end
