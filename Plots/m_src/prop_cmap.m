function cm = prop_cmap2(n0,n1)
  
  if nargin == 0
    n0 = 128;
    n1 = 128;
  end
 
 
  o0 = ones(n0,1);
  c0 = linspace(0,.5, n0)';
 
  o1 = ones(n1,1);
  c1 = linspace(.5,0, n1)';
 
  cm =   [ o0 c0 c0 ; 
           c1 o1 c1];
 
  colormap(cm);
  