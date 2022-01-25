function new_x = cover_complete(x, range, delta, delta_is_float)
%%
% new_x = cover_complete(x, range, delta, delta_is_int)
%
% new_x contains one and only one element of x in grd.
%
% delta_is_float is false by default, but is used only to enforce float
% in case delta is int and should be float. 
% 
% i.e. if delta is a float and delta_is_float is false, it will be ignored
% and delta will not be changed. 
%


if ~exist('range','var')||isempty(range)
    range = [min(x), max(x)];
end

if ~exist('delta', 'var')
    delta = 20; 
end

if ~exist('delta_is_float', 'var')
    delta_is_float = false;
end

if mod(delta,1)~=0  
    delta_is_float = true;
end

if delta_is_float
    grd = range(1):delta:range(2);
else
    grd = linspace(range(1), range(2), delta);
end

x = [x grd]; % augment x with grid

Glow = repmat(grd', 1, size(x,2));
Gup = Glow(2:end,:);
Gup(end+1,:) = inf;

X = repmat(x, size(grd,2),1);
Xin = (X>Glow)&(X<Gup);

new_x = grd;

for i_grd = 1:numel(grd)
    ifound = find(Xin(i_grd,:), 1);  % find the first occurence of x in cell i_grd    
    if ~isempty(ifound)
        new_x(end+1) = x(ifound);
    end
end


