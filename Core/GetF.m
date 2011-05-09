function f = GetF(Sys,t,pts)

  InitSystem(Sys);  
  f = cvm(50,t,pts);