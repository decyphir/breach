function h = plot_histo(M,S,iX, props, iP)
  
  h = figure;    
  nb_histo = numel(iX)+numel(props);
  
  % y labels
  ytick_labels = {};
  for k = 1:numel(iP)
    ylabel = S.ParamList{iP(k)};
    if (iP(k)<= S.DimX)
      ylabel = [ylabel '(0)'];
    end
    ytick_labels = {ytick_labels{:}, ylabel };                    
  end
  
  
  % plotting sensitivities of variables
  
  nh = min(3,nb_histo);
  
  for i = 1:numel(iX)
    subplot(ceil(nb_histo/3),nh,i);
    barh(M(i,:));          
    set(gca, 'YTick', 1:numel(iP), 'YTickLabel',  ytick_labels);    
    axis tight;
    grid on;
    hy = get(gca, 'ylabel');
    set(hy, 'Interpreter','none');        
    st = ['S(' S.ParamList{i} '[t])'];
    title(st, 'Interpreter','none');
    
  end
  
  % plotting sensitivities of properties TODO 
  for i = numel(iX)+1:nb_histo
    subplot(ceil(nb_histo/3),nh,i+numel(iX));
    barh(M(i,:));          
    set(gca, 'YTick', 1:numel(iP), 'YTickLabel',  ytick_labels);    
    hy = get(gca, 'ylabel')
    set(hy, 'Interpreter','none');        
  end
  
  fig_resize(gcf, nh,ceil(nb_histo/3));

  