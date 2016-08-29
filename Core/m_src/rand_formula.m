function [ st_breach, st_taliro, preds ] = rand_formula(n, op, preds)
%RAND_FORMULA Creates a random formula of size n, minimizing the height of the parse tree  
%
% Synopsis: [ st_breach, st_taliro, preds ] = rand_formula(n, op)
% 
% by default op is a subset  {'and', 'ev', 'alw', 'until', 'evI', 'alwI',
% 'untilI'} (contains all operators if not specified)
% 
% returns the random formula with predicates mu1, ..., mun in both Breach and
% TaLiRo formats. 
% 
% SEE ALSO rand_formula_peigne 

  if nargin == 1
    op = {'or', 'and', 'ev', 'alw', 'until', 'evI', 'alwI', 'untilI'};
    preds = {'mu'};
  end
  
  if nargin == 2
    preds = {'mu'};
  end
 
  if (n==0) 
    st_breach = genvarname('mu',{'mu',  preds{:}});
    st_taliro = st_breach;
    preds = {'mu',preds{:}, st_breach};
    return;
  else
      opi = op{randi(numel(op))};
  end
  
  I = [0 0];
  while I(1)==I(2)
      I = sort(20*rand([1 2]));
  end
  
  switch opi
    case 'and'
      [st_breach1, st_taliro1, preds] = rand_formula(floor((n-1)/2), op, preds);
      [st_breach2, st_taliro2, preds] = rand_formula(ceil((n-1)/2), op, preds);
      st_breach  =  ['(' st_breach1  ') and (' st_breach2 ')'];
      st_taliro  =  ['(' st_taliro1  ') /\ (' st_taliro2 ')'];

    case 'or'
      [st_breach1, st_taliro1, preds] = rand_formula(floor((n-1)/2), op, preds);
      [st_breach2, st_taliro2, preds] = rand_formula(ceil((n-1)/2), op, preds);
      st_breach  =  ['(' st_breach1  ') or (' st_breach2 ')'];
      st_taliro  =  ['(' st_taliro1  ') \/ (' st_taliro2 ')'];
   
   case 'ev' 
    [st_breach1, st_taliro1, preds] = rand_formula(n-1, op, preds);
    st_breach  =  ['ev (' st_breach1 ')']; 
    st_taliro  =  ['<> (' st_taliro1 ')']; 

   case 'alw' 
    [st_breach1, st_taliro1, preds] = rand_formula(n-1, op, preds);
    st_breach  =  ['alw (' st_breach1 ')']; 
    st_taliro  =  ['[] (' st_taliro1 ')']; 

   case 'evI' 
    [st_breach1, st_taliro1, preds] = rand_formula(n-1, op, preds);
    st_breach  =  ['ev_[' num2str(I(1),3) ',' num2str(I(2),3)   ']  (' st_breach1 ')']; 
    st_taliro  =  ['<>_[' num2str(I(1),3) ',' num2str(I(2),3)  ']  (' st_taliro1 ')']; 
    
   case 'alwI' 
    [st_breach1, st_taliro1, preds] = rand_formula(n-1, op, preds);
    st_breach  =  ['alw_[' num2str(I(1),3) ',' num2str(I(2),3)  ']  (' st_breach1 ')']; 
    st_taliro  =  ['[]_[' num2str(I(1),3) ',' num2str(I(2),3)  ']  (' st_taliro1 ')']; 
 
   case 'until'
    [st_breach1, st_taliro1, preds] = rand_formula(floor((n-1)/2), op, preds);
    [st_breach2, st_taliro2, preds] = rand_formula(ceil((n-1)/2), op, preds);
    st_breach  =  ['(' st_breach1  ') until (' st_breach2 ')'];  
    st_taliro  =  ['(' st_taliro1  ') U (' st_taliro2 ')'];    
    
   case 'untilI'
    [st_breach1, st_taliro1, preds] = rand_formula(floor((n-1)/2), op, preds);
    [st_breach2, st_taliro2, preds] = rand_formula(ceil((n-1)/2), op, preds);
    st_breach  =  ['(' st_breach1  ') until_[' num2str(I(1),3) ',' num2str(I(2),3)  '] (' st_breach2 ')'];  
    st_taliro  =  ['(' st_taliro1  ') U_[' num2str(I(1),3) ',' num2str(I(2),3)  '] (' st_taliro2 ')'];    

  end
  if (strcmp(preds{1},'mu'))
    preds = {preds{2:end}};    
  end
end
