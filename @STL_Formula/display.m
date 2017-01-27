function st = display(phis)
%DISPLAY displays a set of formulas
% 
% Synopsis: st = display(phis)
% 
% Input:
%  - phis : an array of STL formula
%

for ii = 1:numel(phis)
    st = disp(phis(ii),1);
    fprintf('\n');
    fprintf(st);
    fprintf('\n');
end

fprintf('\n');

end
