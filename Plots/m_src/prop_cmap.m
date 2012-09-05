function cm = prop_cmap(varargin) 

  switch nargin 
   case 0
    n0 = 128;
    n1 = 128;
   case 1
    val= varargin{1};
    M = abs(max(val));
    m = abs(min(val));
    n0= ceil(m/(m+M)*256);
    n1= ceil(M/(m+M)*256);
   case 2
    n0 = varargin{1};
    n1 = varargin{2};          
  end
   
  r0 = linspace(1,0.5, n0)'; %ones(n0,1);
  v0 = linspace(0,.5, n0)';
  b0 = zeros(n0, 1);
  
  r1 = linspace(0.5,0, n1)'; % o1 = ones(n1,1);
  v1 = linspace(.5,1, n1)'; 
  b1 = zeros(n1, 1);
  
  cm =   [ r0 v0 b0 ; 
           r1 v1 b1];
 
  colormap(cm);
  