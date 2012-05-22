function f = GetF(Sys,t,pts)
% GETF returns the value of rhs function of the ODE dynamics
%  
% Synopsis:   f =  GetF(Sys, t, pts) 
%  
%    where Sys is a system, t a time value (can be an array), pts 
%
%    This function can be useful for debugging a new dynamics. 
%
%    Example:    
%    
%     CreateSystem;
%     f = GetF(Sys, 0, Sys.p) % tests the rhs function for nominal parameters of Sys
%  
  
  InitSystem(Sys);  
  f = cvm(50,t,pts);