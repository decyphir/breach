function cm = prop_cmap(varargin) 

  nb_col = 256;
  
  switch nargin 
   case 0
    n0 = nb_col/2;
    n1 = nb_col/2;
    cm = bicolormap(n0,n1);
   case 1
    val= varargin{1};
    
    M = max(val);
    m = min(val);
        
    if (M<0) % all values negative
      
      r = linspace(0.5, 1, nb_col)';
      v = zeros(nb_col,1);
      b = zeros(nb_col,1);            
      
      cm = [r v b];
    elseif (m>0) % all values positive
      
      r = zeros(nb_col,1);
      v = linspace(0.5, 1, nb_col)';
      b = zeros(nb_col,1);           
      cm = [r v b ];
    else % normal situation
      
      M = abs(max(val));
      m = abs(min(val));
      
      n0= ceil(m/(m+M)*nb_col);
      n1= ceil(M/(m+M)*nb_col);
      
      cm = bicolormap(n0,n1);
    end
   case 2
    n0 = varargin{1};
    n1 = varargin{2};      
    cm = bicolormap(n0,n1);
    
  end
   
  colormap(cm);

function cm = bicolormap(n0,n1)
  
  r0 = linspace(1,0.5, n0)'; %ones(n0,1);
  v0 = linspace(0,.5, n0)';
  b0 = zeros(n0, 1);
  
  r1 = linspace(0.5,0, n1)'; % o1 = ones(n1,1);
  v1 = linspace(.5,1, n1)'; 
  b1 = zeros(n1, 1);
  
  cm =   [ r0 v0 b0 ; 
           r1 v1 b1];
 
