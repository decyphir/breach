function  disp_cover_opts(opts)

st = '';
if isfield(opts,'params')
    params = fieldnames(opts.params);
    st = sprintf('Parameters:\n-------------\n');
    for ip= 1:numel(params)
        param_ip= params{ip};
        range_ip = opts.params.(param_ip).range;
        grid_ip = opts.params.(param_ip).grid;
        st_ip = [param_ip ' in [' num2str(range_ip(1)) ', ' num2str(range_ip(2)) '] with ' num2str(numel(grid_ip)) ' bins.'];
        st = sprintf([st '%s\n'], st_ip);
    end
end

if isfield(opts, 'projections')
    st_proj = sprintf('\nProjections:\n-------------\n');
    num_dim  = numel(opts.projections);
    for dim  = 1:num_dim
        st_dim=sprintf([num2str(dim) 'd: ']);
        pref= repmat(' ', 1, numel(st_dim)+3);
        this_dim_proj = opts.projections{dim};
        % first projection
        this_proj = this_dim_proj{1};
        this_proj_st = this_proj{1};
        for j = 2:dim
            this_proj_st = sprintf([this_proj_st '__x__%s'],this_proj{j});
        end
        st_dim  = append_or_linebreak(st_dim, this_proj_st,pref);
        
        % rest of it
        for i = 2:numel(opts.projections{dim})
            this_proj = this_dim_proj{i};
            this_proj_st = this_proj{1};
            for j = 2:dim
                this_proj_st = sprintf([this_proj_st '__x__%s'],this_proj{j});
            end
            st_dim  = append_or_linebreak([st_dim ', '], this_proj_st,pref);
        end
        st_proj = sprintf([st_proj st_dim '\n']);
    end
    st = [st st_proj];
end

if isfield(opts,'signals')
    signals = fieldnames(opts.signals);
    st_sigs = sprintf('\nSignals:\n-------\n');
    for ip= 1:numel(signals)
        signal_ip= signals{ip};
        range_ip = opts.signals.(signal_ip).range;
        grid_ip = opts.signals.(signal_ip).grid;
        st_ip = [signal_ip ' in [' num2str(range_ip(1)) ', ' num2str(range_ip(2)) '] with ' num2str(numel(grid_ip)) ' bins.'];
        st_sigs = sprintf([st_sigs '%s\n'], st_ip);
    end
    st  = [st st_sigs];
end

if isfield(opts, 'signals_freqs')
   signals_freqs = fieldnames(opts.signals_freqs);
   if ~exist('st_sigs','var')
       
       st_sigs = sprintf('\nSignals:\n-------\n');
   end
    for ip= 1:numel(signals_freqs)
        signal_ip= signals_freqs{ip};
        range_ip = opts.signals_freqs.(signal_ip).range;
        grid_ip = opts.signals_freqs.(signal_ip).grid;
        st_ip = [signal_ip ' in [' num2str(range_ip(1)) ', ' num2str(range_ip(2)) '] with ' num2str(numel(grid_ip)) ' bins.'];
        st_sigs = sprintf([st_sigs '%s\n'], st_ip);
    end
    st  = [st st_sigs];
end

if isfield(opts, 'signals_projections')
    st_proj = sprintf('\nSignals Projections:\n-------------\n');
    num_dim  = numel(opts.signals_projections);
    for dim  = 1:num_dim
        st_dim=sprintf([num2str(dim) 'd: ']);
        pref= repmat(' ', 1, numel(st_dim)+3);
        this_dim_proj = opts.signals_projections{dim};
        % first projection
        this_proj = this_dim_proj{1};
        this_proj_st = this_proj{1};
        for j = 2:dim
            this_proj_st = sprintf([this_proj_st '__x__%s'],this_proj{j});
        end
        st_dim  = append_or_linebreak(st_dim, this_proj_st,pref);
        
        % rest of it
        for i = 2:numel(opts.signals_projections{dim})
            this_proj = this_dim_proj{i};
            this_proj_st = this_proj{1};
            for j = 2:dim
                this_proj_st = sprintf([this_proj_st '__x__%s'],this_proj{j});
            end
            st_dim  = append_or_linebreak([st_dim ', '], this_proj_st,pref);
        end
        st_proj = sprintf([st_proj st_dim '\n']);
    end
st  = [st st_proj];    
end


if isfield(opts, 'time_signals_projections')
    st_proj = sprintf('\nTime Signals Projections:\n--------------------------\n');
    num_dim  = numel(opts.time_signals_projections);
    for dim  = 2:num_dim+1        
        st_dim=sprintf([num2str(dim-1) 'd: ']);        
        pref= repmat(' ', 1, numel(st_dim)+3);
        this_dim_proj = opts.time_signals_projections{dim-1};
        % first projection
        this_proj = this_dim_proj{1};
        this_proj_st =  this_proj{1};
        for j = 2:dim
            this_proj_st = sprintf([this_proj_st '__x__%s'],this_proj{j});
        end
        st_dim  = append_or_linebreak(st_dim, this_proj_st,pref);
        
        % rest of it
        for i = 2:numel(opts.time_signals_projections{dim-1})
            this_proj = this_dim_proj{i};
            this_proj_st = this_proj{1};
            for j = 2:dim
                this_proj_st = sprintf([this_proj_st '__x__%s'],this_proj{j});
            end
            st_dim  = append_or_linebreak([st_dim ', '], this_proj_st,pref);
        end
        st_proj = sprintf([st_proj st_dim '\n']);
    end
st  = [st st_proj];    
end

if isfield(opts, 'signals_freqs_projections')
    st_proj = sprintf('\nSignals Frequencies Projections:\n-----------------------------\n');
    num_dim  = numel(opts.signals_freqs_projections);
    for dim  = 1:num_dim        
        st_dim=sprintf([num2str(dim) 'd: ']);        
        pref= repmat(' ', 1, numel(st_dim)+3);
        this_dim_proj = opts.signals_freqs_projections{dim};
        % first projection
        this_proj = this_dim_proj{1};
        this_proj_st =  this_proj{1};
        for j = 2:dim
            this_proj_st = sprintf([this_proj_st '__x__%s'],this_proj{j});
        end
        st_dim  = append_or_linebreak(st_dim, this_proj_st,pref);
        
        % rest of it
        for i = 2:numel(opts.signals_freqs_projections{dim})
            this_proj = this_dim_proj{i};
            this_proj_st = this_proj{1};
            for j = 2:dim
                this_proj_st = sprintf([this_proj_st '__x__%s'],this_proj{j});
            end
            st_dim  = append_or_linebreak([st_dim ', '], this_proj_st,pref);
        end
        st_proj = sprintf([st_proj st_dim '\n']);
    end
st  = [st st_proj];    
end


disp(st);
end