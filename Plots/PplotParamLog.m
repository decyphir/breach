function PplotParamLog(P)
%PPLOTPARAMLOG plots the value of each parameter separately.
% 
% Synopsis: PplotParamLog(P)
% 
% Input:
%  - P : the parameter set for which parameter values are plotted
% 
% Output:
%  - None, but a nice figure
%

figure(gcf);
nParam = size(P.pts,2);
toRemove = floor(nParam*.125); % colormap part to remove to avoid dark colors
my_map = jet(nParam+2*toRemove+1); %size(my_map,2) must be nParam+1
my_map = my_map(1+toRemove:end-toRemove,:);
my_map(1,:) = [1,1,1]; % lowest color is white
colormap(my_map); % define the color map associated to the figure


for ii = 1:numel(P.ParamList)
    param = P.ParamList{ii};
    val_log = log10(GetParam(P,param)); % param value on log scale
    vl_min = min(val_log);
    vl_max = max(val_log);
    if(abs(vl_max-vl_min)<0.1) % in case of a singular point
        if(vl_min>0)
            vl_min = vl_min/1.02;
        else
            vl_min = vl_min*1.02;
        end
        if(vl_max>0)
            vl_max = vl_max*1.02;
        else
            vl_max = vl_max/1.02;
        end
    end
    
    % compute repartition
    N = 15; % == number of spliting intervals
    edge = [-Inf,vl_min+(1:N-1).*((vl_max-vl_min)/N),Inf]; % numel(edge)==N+1
    c = histc(val_log,edge); % c(end) is numel(find(val_log==Inf))==0. Keep it so the plot looks nicer
    
    % define bar x-values
    edge(1) = vl_min;
    edge(end) = vl_max;
    edge = (edge(1:end-1)+edge(2:end))/2; % middle of edge (numel(edge)=N)
    
    % plot
    subplot(numel(P.ParamList),1,ii);
    h = bar(edge,ones(1,N),1);
    delete(findobj('marker','*')) % delete all (unwanted) stars on x axis (maybe not mandatory)
    set(gca,'XLim',[vl_min,vl_max],'YTick',[],'YLim',[0,1]); % so the bar fill all vertical space
    set(gca,'TickDir','out','TickLength',get(gca,'TickLength').*2./3);
    
    % define the colors
    caxis([0,nParam]); % define the range of color value for the current axis
    bar_child = get(h,'Children');
    set(bar_child,'CData',c(1:end-1)) % define the color for each bar
    
    % print ticks
    XTick = ceil(vl_min):floor(vl_max);
    if(10^vl_min<0.96*10^ceil(vl_min)) % if far enougth from last tick
        XTick = [vl_min,XTick];
    end
    if(10^vl_max>1.04*10^floor(vl_max)) % if far enougth from last tick
        XTick = [XTick,vl_max];
    end
    XTickLabel = cellstr(num2str((10.^XTick)','%.1e'))'; % convert to cell array with strings
    % we may improve XTickLabel by showing a "text" based on LaTeX such
    % that if XTick is integer XTickLabel is 10^XTick, and a.10^b otherwise
    set(gca,'XTick',XTick);
    set(gca,'XTicklabel',XTickLabel);
    box off;
    
    ylabel(param,'Interpreter','none','Rotation',0,'HorizontalAlignment','Right');
    
end

% print a colorbar in an other windows.
% All this part is based on matlab default values
figure();
axis off
pos = get(gcf,'Position');
set(gcf,'Position', [pos(1), pos(2), pos(3)/5, pos(4)]);
caxis([0,100]); % define lowest and highest color
colormap(my_map);
%title('Legend');
c = colorbar();
x1=get(gca,'position'); % set colorbar width (code get from matlab website)
x=get(c,'Position');
x(1)=x(1)/2;
x(3)=x(3)*5;
set(c,'Position',x)
set(gca,'position',x1)

end


function PplotParamLog_old(P)
%PPLOTPARAMLOG plots a each parameter separately. For each parameter, a
% star is drawn at the value of the parameter
%

for ii = 1:numel(P.ParamList)
    param = P.ParamList{ii};
    val = GetParam(P,param);

    subplot(numel(P.ParamList),1,ii);
    plot(val,ones(1,numel(val)),'*');
    set(gca,'XScale','log','YTick',[]);
    ylabel(P.ParamList{ii},'Interpreter','none','Rotation',0,'HorizontalAlignment','Right');
    
end


end
