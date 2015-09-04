function J = GetJf(Sys,t,pts)

  InitSystem(Sys);  
  J = cvm(51,t,pts);
