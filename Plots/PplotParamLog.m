function [params, h] = PplotParamLog(P, varargin)
%PPLOTPARAMLOG_COLOR plots the value of each parameter separately in the
% current figure.
% 
% Synopsis: [params, h] = PplotParamLog(P[, params[, Pall[, method[, title[, common_scale]]]]])
% 
% Inputs:
%  - P            : the parameter set for which parameter values are
%                   plotted
%  - params       : (Optional, default or empty=P.ParamList) names or
%                   indexes in P of the parameter to be plotted. Unvalid
%                   parameter names or indexes are not plotted
%  - Pall         : (Optional, default or empty=min and max of values)
%                   parameter set used to define the lower and upper bound
%                   of the plot of each parameter. It must contain only one
%                   parameter vector. It is used only for parameter which
%                   are uncertain parameter of Pall. Otherwise, the minimal
%                   and maximal value of parameter in P is considered.
%  - method :       (Optional, default or empty='color') string describing
%                   the method used to plot parameter values. It must be
%                   either 'color' or 'star'.
%  - title        : (Optional, default='') string describing the title
%                   shown at the top of the figure. It is interpreted with
%                   TeX interpretor.
%  - common_scale : (Optional, default=false) boolean true if all plot must
%                   span over the same range.
% 
% Outputs:
%  - params : a cell array containing the name of plotted parameter
%  - h      : array of handle containing handle of each figure shown in the
%             order of appearance. If using 'star' method, it contains the
%             handle of the only shown figure. If using 'color' method, it
%             contains in the first position, the handle of the figure, and
%             in second position, the handle of the legend.
% 
% Example (Lorentz84):
%   CreateSystem;
%   P = CreateParamSet(Sys,{'a','b'},[1,10;10,1000]);
%   P = SetParam(P,{'x0','x1','x2'},[1;1;1]);
%   P = Refine(P,10); % 100 parameter vectors generated
%   figure();
%   PplotParamLog(P); % plot repartition (on log scale)
%   
%   P = CreateParamSet(Sys,{'a','b'},[1,10;10,1000]);
%   P = SetParam(P,{'x0','x1','x2'},[1;1;1]);
%   P = RandomLogRefine(P,100); % 100 parameter vectors generated
%   figure();
%   PplotParamLog(P,{'x0','a','b'}); % homogeneous repartition (log scale)
%   
%   P = CreateParamSet(Sys,{'x0','x1','x2'},[1,10;10,300;4,1000]);
%   P2 = Refine(P,20);
%   figure();
%   PplotParamLog(P2,{'x0','x1','x2'},P,'','',true);
%   figure();
%   PplotParamLog(P2,{'x0','x1','x2'},P,'star','',true); % star version
% 
%See also SplotVar SplotTraj
%

if(nargin>=4)
    method = varargin{3};
    varargin(3) = []; % remove the method from varargin
    if strcmpi(method,'star')
        [params,h] = PplotParamLog_star(P,varargin{:});
    else
        [params,h] = PplotParamLog_color(P,varargin{:});
    end
else
    [params,h] = PplotParamLog_color(P,varargin{:});
end

end

function [params, h] = PplotParamLog_color(P, params, Pall, title_str, common_scale)
%PPLOTPARAMLOG_COLOR plots the value of each parameter separately.
% 
% Synopsis: [params, h] = PplotParamLog(P[, params[, Pall[, title_str[, common_scale]]]])
% 
% Inputs:
%  - P            : the parameter set for which parameter values are
%                   plotted
%  - params       : cell array of string or indexes of parameters to plot
%  - Pall         : parameter set describing the interval of the parameter
%  - title_str    : string containing the graph title
%  - common_scale : boolean set to true if all plot are on the same scale
% 
% Outputs:
%  - params : name of plotted parameters
%  - h      : array containing the handle of the figure and the handle of
%             the legend
%

if ~exist('params','var')
    params = {};
elseif ischar(params)
    params = {params};
elseif ~iscell(params)
    params = P.ParamList(params);
end
valid = ismember(params,P.ParamList);
params = params(valid);
if isempty(params)
    params = P.ParamList;
end

if ~exist('common_scale','var')
    common_scale = false;
end

nParam = numel(params); % number of parameters to plot
nParamVect = size(P.pts,2); % number of parameter vectors
if(common_scale)
    % variables used iff common_scale if true
    min_plot = Inf;
    max_plot = -Inf;
    empty_box = []; % plot not to redraw
end

% define color map
h1 = figure(gcf);
toRemove = floor(nParamVect*.125); % colormap part to remove to avoid dark colors
my_map = jet(nParamVect+2*toRemove+1);
my_map = my_map(1+toRemove:end-toRemove,:); %size(my_map,2) must be nParamVect+1
my_map(1,:) = [1,1,1]; % lowest color is white
colormap(my_map); % define the color map associated to the figure

for ii = 1:nParam
    subplot(nParam,1,ii); % one subplot for each parameter
    param = params{ii};
    val_log = log10(GetParam(P,param)); % param value on log scale
    
    % clean val_log values
    vl_real = arrayfun(@(x) isreal(x{:}), num2cell(val_log)); % remove complex
    if any(~vl_real)
        warning('PplotParamLog:negativeValue',...
            'At least one value for parameter %s is negative, it is removed.',param);
        val_log = val_log(vl_real);
    end
    vl_inf = isinf(val_log); % remove inf
    if any(vl_inf)
        warning('PplotParamLog:infiniteValue',...
            'One (or more) value for parameter %s is 0 or Inf, it is removed.',param);
        val_log = val_log(~vl_inf);
    end
    
    % check validity of parameter
    if isempty(val_log)
        box on;
        set(gca,'XTick',[],'YTick',[]); % display empty graph
        ylabel(param,'Interpreter','none','Rotation',0,'HorizontalAlignment','Right');
        if(common_scale)
            empty_box = [empty_box,ii]; % don't redraw this one
        end
        continue;
    end
    
    % define parameter upper and lower bounds
    if( exist('Pall','var') && ~isempty(Pall) && ~isempty(GetEpsi(Pall,param)) )
        ptsAll = GetParam(Pall,param);
        epsiAll = GetEpsi(Pall,param);
        vl_min = log10(ptsAll-epsiAll);
        vl_max = log10(ptsAll+epsiAll);
        if(isinf(vl_min) || ~isreal(vl_min))
            warning('PplotParamLog:badLowerLimit',...
                'The lower bound for parameter %s is <=0 or Inf, set to lowest valid value of param set.',param);
            vl_min = min(val_log);
        end
        if(isinf(vl_max) || ~isreal(vl_max))
            warning('PplotParamLog:badLowerLimit',...
                'The upper bound for parameter %s is <=0 or Inf, set to highest valid value of param set.',param);
            vl_max = max(val_log);
        end
        if(min(val_log)<vl_min)
            warning('PplotParamLog:badLowerBound',...
                'Some parameter values for %s (aka: %g) are lower than the interval (%g).',param,10^min(val_log),10^vl_min);
        end
        if(max(val_log)>vl_max)
            warning('PplotParamLog:badUpperBound',...
                'Some parameter values for %s (aka: %g) are higher than the interval (%g).',param,10^max(val_log),10^vl_max);
        end
    else
        if(exist('Pall','var') && ~isempty(Pall)) % here isempty(GetEpsi(Pall,param)) holds
            warning('PplotParamLog:uncompletePall',...
                'The parameter set Pall does not define the interval for parameter %s. Using valid values to define interval.',param);
        end
        vl_min = min(val_log);
        vl_max = max(val_log);
    end
    if(abs(vl_max-vl_min)<0.1) % in case of a singular point
        vl_min = vl_min - 0.05;
        vl_max = vl_max + 0.05;
    end
    
    % save highest and lowest bounds
    if(common_scale)
        min_plot = min(min_plot,vl_min);
        max_plot = max(max_plot,vl_max);
    end
    
    % compute repartition
    N = 15; % == number of spliting intervals
    edge = [-Inf,vl_min+(1:N-1).*((vl_max-vl_min)/N),Inf]; % numel(edge)==N+1
    histo = histc(val_log,edge); % histo(end) is numel(find(val_log==Inf)) is 0
    
    % define bar x-values
    edge(1) = vl_min;
    edge(end) = vl_max;
    edge = (edge(1:end-1)+edge(2:end))/2; % middle of edge (numel(edge)=N)
    
    % plot
    h = bar(edge,ones(1,N),1);
    delete(findobj('marker','*')) % delete all (unwanted) stars on x axis (maybe not mandatory)
    set(gca,'XLim',[vl_min,vl_max],'YLim',[0,1]); % make the bar fills all vertical space
    
    % define the colors
    caxis([0,nParamVect]); % define the range of color value for the current axis
    bar_child = get(h,'Children');
    set(bar_child,'CData',histo(1:end-1)) % define the color for each bar
    
    % define and print ticks
    set(gca,'YTick',[]);
    set(gca,'TickDir','out','TickLength',get(gca,'TickLength').*2./3);
    tick_step = ceil((floor(vl_max)-ceil(vl_min))/13); % 1 if [2,14] ticks, 2 if [15,27] ticks
    tick_step = max(tick_step,1); % tick_step may be 0
    XTick = ceil(vl_min):tick_step:floor(vl_max);
    if isempty(XTick)
        XTick = [vl_min,vl_max];
    else
        if(ceil(vl_min)-vl_min > 0.04*(vl_max-vl_min)) % if far enougth from last tick
            XTick = [vl_min,XTick];
        end
        if(vl_max-XTick(end) > 0.04*(vl_max-vl_min)) % if far enougth from last tick
            XTick = [XTick,vl_max];
        end
    end
    XTickLabel = cellstr(num2str((10.^XTick)','%.1e'))'; % convert to cell array with strings
    % we may improve XTickLabel by showing a "text" based on LaTeX such
    % that if XTick is integer XTickLabel is 10^XTick, and a.10^b otherwise
    set(gca,'XTick',XTick);
    set(gca,'XTicklabel',XTickLabel);
    box off;
    
    ylabel(param,'Interpreter','none','Rotation',0,'HorizontalAlignment','Right');
    if(ii==1 && exist('title_str','var'))
        title(title_str);
    end
end

% rescale plots
if(common_scale && numel(empty_box)<nParam)
    for ii = 1:nParam
        if(~ismember(ii,empty_box))
            subplot(nParam,1,ii);
            xl = get(gca,'XLim');
            set(gca,'XLim',[min_plot,max_plot]);
            rectangle('Position',[min_plot,0,max_plot-min_plot,1],'FaceColor','none','EdgeColor','k'); % homogeneous border
            if(xl(1)-min_plot>0)
                rectangle('Position',[min_plot,0,xl(1)-min_plot,1],'FaceColor',[.75,.75,.75],'EdgeColor','k'); % gray rectangle
            end
            if(max_plot-xl(2)>0)
                rectangle('Position',[xl(2),0,max_plot-xl(2),1],'FaceColor',[.75,.75,.75],'EdgeColor','k'); % gray rectangle
            end
        end
    end
end

% print a colorbar in an other windows.
% All this part is based on matlab default values
h2 = figure();
axis off
pos = get(gcf,'Position');
set(gcf,'Position', [pos(1), pos(2), pos(3)/5, pos(4)]);
caxis([0,100]); % define lowest and highest color
colormap(my_map);
%title('Legend');
cb = colorbar();
x1=get(gca,'position'); % set colorbar width (code get from matlab website)
x=get(cb,'Position');
x(1)=x(1)/2;
x(3)=x(3)*5;
set(cb,'Position',x)
set(gca,'position',x1)

h = [h1,h2]; % array of handle

end


function [params,h] = PplotParamLog_star(P, params, Pall, title_str, common_scale)
%PPLOTPARAMLOG_STAR plots a each parameter separately. For each parameter,
% a star is drawn at the value of the parameter
% 
% Synopsis: [params, h] = PplotParamLog_star(P[, params[, Pall[, title_str[, common_scale]]]])
%

if ~exist('params','var')
    params = {};
elseif ischar(params)
    params = {params};
elseif ~iscell(params)
    params = P.ParamList(params);
end
valid = ismember(params,P.ParamList);
params = params(valid);
if isempty(params)
    params = P.ParamList;
end

if ~exist('common_scale','var')
    common_scale = false;
end

nParam = numel(params);
if(common_scale)
    % variables used iff common_scale if true
    min_plot = Inf;
    max_plot = -Inf;
end

h = figure(gcf);
for ii = 1:nParam
    param = params{ii};
    val = GetParam(P,param);

    subplot(nParam,1,ii);
    plot(val,ones(1,numel(val)),'*');
    if(exist('Pall','var')&&~isempty(Pall))
        epsiAll = GetEpsi(Pall,param);
        if ~isempty(epsiAll)
            ptsAll = GetParam(Pall,param);
            set(gca,'XLim',[ptsAll-epsiAll,ptsAll+epsiAll]);
        end
    end
    set(gca,'XScale','log','YTick',[]);
    ylabel(param,'Interpreter','none','Rotation',0,'HorizontalAlignment','Right');
    if(ii==1 && exist('title_str','var'))
        title(title_str);
    end
    
    if(common_scale)
        xl = get(gca,'XLim');
        min_plot = min(min_plot,xl(1));
        max_plot = max(max_plot,xl(2));
    end
    
end

if(common_scale)
    for ii = 1:nParam
        subplot(nParam,1,ii);
        xl = get(gca,'XLim');
        yl = get(gca,'YLim');
        set(gca,'XLim',[min_plot,max_plot]);
        if(xl(1)-min_plot>0)
            rectangle('Position',[min_plot,yl(1),xl(1)-min_plot,yl(2)-yl(1)],'FaceColor',[.75,.75,.75],'EdgeColor','none');
        end
        if(max_plot-xl(2)>0)
            rectangle('Position',[xl(2),yl(1),max_plot-xl(2),yl(2)-yl(1)],'FaceColor',[.75,.75,.75],'EdgeColor','none');
        end
    end
end


end
