function p = highlight_interval(ax, interval, color, alpha)

y=get(ax,'YLim');
p= patch('XData',[interval fliplr(interval)],'YData',[y(1) y(1) y(2) y(2)],'FaceColor',color,'FaceAlpha',alpha);
 