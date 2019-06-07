R = MyReq1;

MyReq1.PlotDiagnostic(phi1);


return;
%%
robustness_map = R.robustness_map;
diag_map = R.diag_map;
formula_names_map = R.formula_names_map;
nb_plots = robustness_map.size(1);
keys = robustness_map.keys();


if (verdict)
    color = 'g';
else
    color = 'r';
end

for i=1:nb_plots
    id = keys{i};
    h = F.AddAxes(); 
    axis(h);
    hold on;
    grid on;
    formula_name = formula_names_map(id);
    title(formula_name, 'Interpreter', 'none');
    signal = robustness_map(id);
    stairs(signal.times, signal.values);
    
    ylim = get(h, 'YLim');
    ylim_bot = ylim(1);
    ylim_top = ylim(2);
    
    implicant = diag_map(id);
    size = implicant.getIntervalsSize();
    for j=1:size
        interval = implicant.getInterval(j);
        x = interval.begin;
        y = interval.end;
        if (x == y)
            line([x x],[ylim_bot ylim_top],'Color',color);
        elseif (y > x)
            p = patch([x y y x], [ylim_bot ylim_bot ylim_top ylim_top], color);
            alpha(p, 0.05);
            set(p,'EdgeColor','none');
        end
    end
    
    samples = implicant.getSignificantSamples();
    for j=1:length(samples)
        sample = samples(j);
        hold on;
        plot(sample.time, sample.value, 'x');
    end
end