function clim = sym_clim(val)
  
  M = max(max(val));
  m = min(min(val));
  
  switch (sign(M*m))
   case 1
    if m<0
      clim= [m -m];
    else
      clim = [-M M];
    end
   case 0
     if m==0
       clim = [ -M M];
     else
       clim = [m -m];
     end
   case -1 
    if abs(m)>M
      clim = [m -m];
    else
      clim = [-M M];
    end
  end
  

