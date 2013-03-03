function h = plot_histo(M, S, iX, props, iP)
%PLOT_HISTO plots sensitivity histograms
%
% Synopsis : h = plot_histo(M, S, iX, props, iP)
%
% Inputs:
%  - M     : the matrix containing data
%  - S     : the parameter set
%  - iX    : 
%  - props : 
%  - iP    : indexs of parameters
%

h = figure;
nb_histo = numel(iX)+numel(props);

% y labels
ytick_labels = cell(1,numel(iP));
for ii = 1:numel(iP)
    ylabel = S.ParamList{iP(ii)};
    if (iP(ii)<= S.DimX)
        ylabel = [ylabel '(0)']; %#ok<AGROW>
    end
    ytick_labels(ii) = {ylabel};
end


% plotting sensitivities of variables

nh = min(3,nb_histo);

for ii = 1:numel(iX)
    subplot(ceil(nb_histo/3),nh,ii);
    barh(M(ii,:));
    set(gca, 'YTick', 1:numel(iP), 'YTickLabel', ytick_labels);
    axis tight;
    grid on;
    hy = get(gca, 'ylabel');
    set(hy, 'Interpreter', 'none');
    st = ['S(' S.ParamList{iX(ii)} '[t])'];
    title(st, 'Interpreter', 'none');
    
end

% plotting sensitivities of properties TODO
for ii = numel(iX)+1:nb_histo
    subplot(ceil(nb_histo/3),nh,ii);
    barh(M(ii,:));
    set(gca, 'YTick', 1:numel(iP), 'YTickLabel', ytick_labels);
    hy = get(gca, 'ylabel');
    set(hy, 'Interpreter', 'none');
    %TODO define the title
end

fig_resize(gcf, nh, ceil(nb_histo/3));

end
