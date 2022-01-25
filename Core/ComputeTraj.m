function Pf = ComputeTraj(Sys, P0, tspan, u)
%COMPUTETRAJ computes trajectories for a system given initial conditions
% and parameters
%
% Synopsis:   Pf = ComputeTraj(Sys,P0,[tspan ,u])
%
% Inputs:
% -  Sys   : System (needs to be compiled)
% -  P0    : Initial conditions and params given in a parameter set or in
%            an array of size Sys.DimP x size(P0,2). size(P0,2) is the
%            number of trajectories that will be computed.
% -  tspan : Interval of the form [t0, tf]. Fixed time instants can also
%            be specified tspan = [t0 t1 ... tN]; if absent, uses
%            Sys.tspan.
% -  u     : is a structure array of nb_traj inputs. u(i) must be a struct
%            with fields
%             - params_idx : indicates which parameters are made time
%                 dependant
%             - tin : times when the input changes
%             - values : values of the parameters
%
% Output:
%  -  Pf   Parameter set augmented with the field traj containing
%          computed trajectories if the input is a param set. The field
%          traj_ref is filled. If the P0 is an  array of parameter values,
%          then Pf is an array of trajectories.
%
% Examples (Lorentz84):
%   CreateSystem;
%   P = CreateParamSet(Sys,'a',[1,2]);
%
%   P1 = Refine(P,2);
%   P1 = ComputeTraj(Sys,P1,0:0.1:10);
%   P1 = ComputeTraj(Sys,P1,0:0.1:10); % Here, nothing shows because nothing happens
%
%   P2 = SetParam(P,'paramProp',2);
%   P2 = SAddUncertainParam(P2,'paramProp');
%   P2 = Refine(P2,2);
%   P2 = ComputeTraj(Sys,P2,0:0.1:10);
%   P2.traj_ref  % should be [1 2 1 2]
%
%See also CreateParamSet Sselect SConcat SPurge
%

if isfield(Sys,'Verbose')
    Verbose = Sys.Verbose;
else
    Verbose=1;
end

if nargin==2
    tspan=Sys.tspan;
end

% checks if we have a parameter set or an array of parameter values

output_trajs = 0; % is 1 if the output must be the array of trajectories

if ~isstruct(P0)
    % We create a parameter set
    if(size(P0,1) ~= Sys.DimP)
        if(size(P0,2) == Sys.DimP) % be smart, try transpose in case it works
            P0 = P0';
        else
            error('ComputTraj:S0DimensionError',...
                'Second argument must be a parameter set or be of dimension Sys.DimP x nb_traj.')
        end
    end
    output_trajs = 1;
    
    pts = P0;
    P0 = CreateParamSet(Sys,1);
    P0.pts = pts;
    P0.epsi = ones(1, size(pts,2));
    P0.traj_to_compute = 1:size(pts,2);
    P0.traj_ref = zeros(1,size(pts,2));
    
end

% checks for an initialization function
if isfield(Sys, 'init_fun')
    P0 = Sys.init_fun(P0);
end

% if no trajectories to compute, we return the param set itself
if(isfield(P0, 'traj_to_compute') && isempty(P0.traj_to_compute))
    if isfield(P0, 'traj')&&isfield(P0, 'traj_ref') 
        for ipts = 1:size(P0.pts,2)
            if ~isequal( tspan(end), P0.traj{P0.traj_ref(ipts)}.time(1,end))
                P0 = rmfield(P0,'traj');
                P0 = Preset_traj_ref(P0);
                break;
            end
        end
        if isempty(P0.traj_to_compute)
            Pf = P0;
            return;
        end
    end
end

if ~isfield(Sys, 'type')
    Sys.type = '';
end

if strcmp(Sys.type,'traces') % No model
    % If system type is only traces, check consistency of params and pts
    Pf = P0;
    for ii = 1:numel(Pf.traj)
        Pf.traj{ii}.param = Pf.pts(1:Pf.DimP,ii)';
    end
 elseif(isfield(P0,'traj_to_compute') &&...
         ~isempty(P0.traj_to_compute) && ~isequal(P0.traj_to_compute,1:size(P0.pts,2))&&... % some traces have already be computed
         isfield(P0, 'traj')&&~isempty(P0.traj)&&isequal(P0.traj{1}.time, tspan))     %  some traces have been computed on the same tspan
     % Here, we assume:
     % 1/ that the index of a param vector is not in traj_to_compute if
     % there is a valid simulation for this param vector
     % 2/ that the field traj_to_compute is ordered and contains unique
     % parameter vectors wrt to system parameter
     Ptmp = Sselect(P0, P0.traj_to_compute);
     Ptmp = ComputeTraj(Sys, Ptmp, tspan);
     
     Pf = P0;
     numTrajP0 = 0;
     if isfield(P0,'traj')
         numTrajP0 = numel(P0.traj); % in case there is already some traj in P0
     end
     [~,~,i_P0] = unique([Ptmp.pts(1:P0.DimP,:),P0.pts(1:P0.DimP,:)]','rows','stable'); % Ptmp(1:P0.DimP,:) are all unique
     for ii = 1:numel(P0.traj_to_compute) % for each newly computed traj
         Pf.traj{numTrajP0+ii} = Ptmp.traj{Ptmp.traj_ref(ii)}; % add it to Pf
         Pf.Xf(1:Pf.DimX,numTrajP0+ii) = Ptmp.Xf(1:Ptmp.DimX,Ptmp.traj_ref(ii));
         i_traj_ref = find(i_P0==ii); % look for indexes of param vector in Pf corresponding to this traj
         i_traj_ref = i_traj_ref(i_traj_ref>numel(P0.traj_to_compute)) - numel(P0.traj_to_compute); % The first ones are Ptmp index, skip them
         Pf.traj_ref(i_traj_ref) = numTrajP0+ii;
     end
     Pf.traj_to_compute = [];
     if(isfield(Sys,'time_mult') && ~isfield(Pf,'time_mult'))
         Pf.time_mult = Sys.time_mult;
     end
     
     return;
end

% From now, we only got unique system-parameter vectors
Pf = P0;
%ipts = 1:size(P0.pts,2);
ipts = P0.traj_to_compute; % all hell let loose, trusting
                           % traj_to_compute and traj_ref... 

ii=0;

switch Sys.type
    case 'Extern'
        model = Sys.name;
        if Verbose==1
            if(numel(ipts)>1)
                rfprintf_reset();
                rfprintf(['Computed ' num2str(ii) '/' num2str(numel(ipts)) ' simulations of ' model])
            end
        end
            
        icount = 0;
        for ii = ipts
            iref = P0.traj_ref(ii);
            if isfield(Sys,'init_u')
                U = Sys.init_u(Sys.ParamList(Sys.DimX-Sys.DimU+1:Sys.DimX), P0.pts(1:Sys.DimP,ii), tspan);
                assignin('base','t__',U.t);
                assignin('base', 'u__',U.u);
            end
            
            [traj.time, traj.X] = Sys.sim(Sys, tspan, P0.pts(:,ii));
            traj.param = P0.pts(1:P0.DimP,ii)';
            
            if isfield(Sys, 'output_gens')
                for io = 1:numel(Sys.output_gens)
                    og = Sys.output_gens{io};
                    % Find in_signals
                    is = FindParam(Sys, og.signals_in);
                    ip = FindParam(P0, og.params);
                    X_in = traj.X(is, :);
                    pts_in = P0.pts(ip,ii);
                    [traj.time, Xout_i] = og.computeSignals(traj.time, X_in, pts_in);
                    traj.X = [traj.X ;Xout_i ];
                end
            end
            
            Pf.traj{iref} = traj;
            Pf.Xf(:,iref) = zeros(Pf.DimX,1);
            icount = icount+1;
                    
            if Verbose==1
                if(numel(ipts)>1)
                    rfprintf(['Computed ' num2str(icount) '/' num2str(numel(ipts)) ' simulations of ' model])
                end
            end
            
        end
        if Verbose==1
            if(numel(ipts)>1)
                fprintf('\n');
            end
        end
        
        Pf.traj_to_compute = [];
        
    case 'Simulink'
        model = Sys.mdl;
        
        if Verbose==1
            if(numel(ipts)>1)
                rfprintf_reset();
                rfprintf(['Computed ' num2str(ii) '/' num2str(numel(ipts)) ' simulations of ' model])
            end
        end
        
        if numel(ipts) == 1
            Verbose=0;
        end
        
        if isfield(Sys, 'Parallel')&&Sys.Parallel&&numel(ipts)>1
            
            %% Parallel computation            
            trajs = cell(1, numel(ipts));            
            batch_size = 100; 
            start_idx = 1;
            end_idx = min(batch_size,numel(ipts));
            num_computed = 0; 
            while start_idx <numel(ipts) 
                
                for idx = start_idx:end_idx
                    f(P0.traj_ref(ipts(idx))) = parfeval(@(ii)task_sim(Sys,P0,tspan,ii), 1, ipts(idx));  % puts ipts(idx) in the queue, such a way that completeIdx will point to the right pos in trajs cell
                end
                
                for idx = start_idx:end_idx
                    % fetchNext blocks until more results are available, and
                    % returns the index into f that is now complete, as well
                    % as the value computed by f.
                    [completedIdx, value] = fetchNext(f);   % fetches something in the queue, which is completedIdx
                    trajs{completedIdx} = value;
                    num_computed =num_computed+1;
                    if Verbose >=1
                        if(numel(ipts)>1)
                            rfprintf(['Computed ' num2str(num_computed) '/' num2str(numel(ipts)) ' simulations of ' model])
                        end
                    end
                end
                start_idx = end_idx+1;
                end_idx = min(end_idx+batch_size,numel(ipts));
            end            
            
            
            Pf.traj = trajs;
            for ii=ipts
                iref = P0.traj_ref(ii);
                Pf.Xf(:,iref) = zeros(Pf.DimX,1);
            end
            
            if Verbose>=1
                if(numel(ipts)>1)
                    fprintf('\n');
                end
            end
            
        else
            %% Serial computation
            trajs = cell(1, numel(ipts));
            icount = 0;
            for ii = ipts
                iref = P0.traj_ref(ii);
                trajs{iref} = task_sim(Sys, P0, tspan, ii);
                if Verbose >=1
                    if(numel(ipts)>1)
                        icount = icount+1;
                        rfprintf(['Computed ' num2str(icount) '/' num2str(numel(ipts)) ' simulations of ' model])
                    end
                end
            end
            Pf.traj = trajs;
            for ii=ipts
                iref = P0.traj_ref(ii);
                Pf.Xf(:,iref) = zeros(Pf.DimX,1);
            end
            if Verbose>=1
                if(numel(ipts)>1)
                    fprintf('\n');
                end
            end
        end
        
        
        Pf.traj_to_compute = [];
        
    otherwise
        
        InitSystem(Sys);
        
        if iscell(tspan)
            if(numel(tspan)==2)
                T = [tspan{1} tspan{2} tspan{2}];  % not really nice.. should be
                                                   % changed some day
            else
                T = cell2mat(tspan);
            end
        else
            T = tspan;
        end
        
        if exist('u','var')
            err = check_u(u);
            if(err~=0)
                error('ComputTraj:ErrorWithU',err);
            end
            
            %This is quite ugly...
            Pf = P0;
            Pf.pts = P0.pts(P0.traj_to_compute,:);
            Pf = cvm(61, Pf, T, u);
            Pf.pts = P0.pts;
            
            Pf.u = u;
            
        else
            Pf = P0;
            Pf.pts = P0.pts(P0.traj_to_compute,:);
            Pf = cvm(61, Pf, T); % <- NM: I would love to know how it works inside!
            Pf.pts = P0.pts;
        end
        
        CVodeFree();
        
        % all trajectories have been computed
        Pf.traj_to_compute = [];
        %Pf.traj_ref = 1:size(P0.pts,2);
        
        if(output_trajs)
            Pf = Pf.traj;
        end
        
end

if(isfield(Sys,'time_mult') && ~isfield(Pf,'time_mult'))
    Pf.time_mult = Sys.time_mult;
end

end

function traj = task_sim(Sys, P0,tspan,ii)

use_caching = isfield(Sys,'DiskCachingFolder')&&(~isempty(Sys.DiskCachingFolder));
do_compute = 1;
p = P0.pts(1:Sys.DimP,ii);
if  use_caching
    hash_traj = DataHash({Sys.ParamList, p, tspan});
    cache_traj_filename = [Sys.DiskCachingFolder filesep 'traj_' hash_traj '.mat'];
    if exist(cache_traj_filename, 'file')
        do_compute = 0;
    else
        do_compute = 1;
    end
    traj = matfile(cache_traj_filename,  'Writable', true);
end

if do_compute
    
    [traj.time, traj.X,traj.status] = Sys.sim(Sys, tspan, P0.pts(:,ii));
    traj.param = P0.pts(1:P0.DimP,ii)';
    
    % compute outputs  - should we do something when status is not right?
    if isfield(Sys, 'output_gens')
        Xout = [];
        for io = 1:numel(Sys.output_gens)
            og = Sys.output_gens{io};
            % Find in_signals
            is = FindParam(Sys, og.signals_in);
            ip = FindParam(P0, og.params);
            X_in = [traj.X(is, :); Xout] ;
            pts_in = P0.pts(ip,ii);
            [traj.time, Xout_i] = og.computeSignals(traj.time, X_in, pts_in);
            Xout = [Xout ; Xout_i];
        end
        traj.X(end-size(Xout,1)+1:end, :) = Xout;
    end
    
    if use_caching % cache new trace
        traj.Properties.Writable = false;
    end
end
end



function err = check_u(u)

err = 0;

if ~isstruct(u)
    err = 'u has to be a structure';
    return;
end

if ~isfield(u,'params_idx')
    err = 'missing field params_idx';
    return;
end

if ~isfield(u,'time')
    err = 'missing field time';
    return;
end

if ~isfield(u,'values')
    err = 'missing field values';
    return;
end

if numel(u.params_idx) ~= size(u.values, 1)
    err = 'numel(u.params_idx) should be equal to  size(u.values, 1)';
    return;
end

if  numel(u.time) ~= size(u.values, 2)
    err = 'numel(u.time) should be equal to size(u.values, 2)';
    return;
end

end
