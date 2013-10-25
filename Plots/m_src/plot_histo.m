function h = plot_histo(M, P, iX, phis, iP)
%PLOT_HISTO plots sensitivity histograms
% 
% Synopsis : h = plot_histo(M, P, iX, phis, iP)
% 
% Inputs:
%  - M    : the matrix containing data
%  - P    : the parameter set
%  - iX   : 
%  - phis : 
%  - iP   : indexes of parameters
%

h = figure;
nb_histo = numel(iX)+numel(phis);

% y labels
ytick_labels = cell(1,numel(iP));
for ii = 1:numel(iP)
    ylabel = P.ParamList{iP(ii)};
    if (iP(ii)<= P.DimX)
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
    st = ['S(' P.ParamList{iX(ii)} '[t])'];
    title(st, 'Interpreter', 'none');
    
end

% plotting sensitivities of properties TODO
for ii = numel(iX)+1:nb_histo
    subplot(ceil(nb_histo/3),nh,ii);
    barh(M(ii,:));
    set(gca, 'YTick', 1:numel(iP), 'YTickLabel', ytick_labels);
    axis tight;
    grid on;
    hy = get(gca, 'ylabel');
    set(hy, 'Interpreter', 'none');
    phi = struct2cell(phis(ii-numel(iX)));
    st = ['S(' phi{1} ')'];
    title(st, 'Interpreter', 'none');
end

fig_resize(gcf, nh, ceil(nb_histo/3));

end

%NM: I copy these two functions from SplotSensiBar
function h = plot_histo3d(Mend, S, iX, props, iP)

figure;
h = bar3(Mend,0.5,'detached');

% Ticks labels

ytick_labels = cell(1,numel(iX)+numel(props));
for ii = 1:numel(iX)
    ytick_labels(ii) = S.ParamList(iX(ii));
end
for ii = 1:numel(props)
    ytick_labels(numel(iX)+ii) = {get_id(props{ii})};
end

xtick_labels = cell(1,numel(iP));
for ii = 1:numel(iP)
    xlabel = S.ParamList{iP(ii)};
    if(iP(ii) <= S.DimX)
        xlabel = [xlabel '(0)']; %#ok<AGROW>
    end
    xtick_labels(ii) = {xlabel};
end

set(gca, 'XTick', 1:size(Mend,2));
set(gca, 'YTick', 1:size(Mend,1));
set(gca, 'XTickLabel',  xtick_labels );
set(gca, 'YTickLabel',  ytick_labels );
axis([0 size(Mend,2)+1 0 size(Mend,1)+1]);

hx = get(gca, 'xlabel');
set(hx, 'Interpreter','none');
hy = get(gca, 'ylabel');
set(hy, 'Interpreter','none');

shading interp;
colormap cool;
%  colorbar;
for i = 1:length(h)
    zdata = get(h(i),'Zdata');
    set(h(i),'Cdata',zdata);
    set(h,'EdgeColor','k');
end

end

function h = plot_histo1(M, S, iX, props, iP)

h = figure;
%nb_histo = numel(iX)+numel(props);

% y labels

ytick_labels = cell(1,numel(iP));
for ii = 1:numel(iP)
    ylabel = S.ParamList{iP(ii)};
    if (iP(ii)<= S.DimX)
        ylabel = [ylabel '(0)'];
    end
    ytick_labels(ii) = {ylabel};
end

% plotting sensitivities of variables


barh(M);
set(gca, 'YTick', 1:numel(iP), 'YTickLabel',  ytick_labels);
hy = get(gca, 'ylabel');
set(hy, 'Interpreter','none');

legend(S.ParamList{iX});

end
