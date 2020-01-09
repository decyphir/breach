function varargout = display(phis)
%DISPLAY displays a set of formulas
% 
% Synopsis: st = display(phis)
% 
% Input:
%  - phis : an array of STL formula
%

varargout = {};
for ii = 1:numel(phis)
    st = disp(phis(ii),1);
    fprintf('\n');
    fprintf(st);
    fprintf('\n');
end

if nargout == 0
    fprintf('\n');
else
    varargout{1} = st;
end

end
