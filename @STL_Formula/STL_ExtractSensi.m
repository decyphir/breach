function params = STL_ExtractSensi(phi)
  
% STL_ExtractSensi find params needing the computation of sensitivities
%

  st = disp(phi);
  [start_idx, end_idx, extents, matches, tokens] = regexp(disp(phi), ['d\{\s*(.+?)\s*\}{(.+?)}\[(.+?)\]']);
  
  params = {};
  for j = 1:numel(matches)
    params = {params{:},tokens{j}{2}};
  end
 
  
