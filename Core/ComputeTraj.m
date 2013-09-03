function Pf = ComputeTraj(Sys, P0, tspan, u)
%  COMPUTETRAJ compute trajectories for a system given initial conditions
%  and parameters
%
%  Synopsis:   Pf = ComputeTraj(Sys,P0,tspan [,u])
%
%   Compute trajectories issued from points in S0 on the
%   time interval tspan
%
%   Inputs:
%
%    -  Sys      System (needs to be compiled)
%
%    -  P0       Initial conditions and params given in a parameter set
%                or in an array of size DimX x nb_traj
%    -  tspan    interval of the form [t0, tf]; Fixed time instants can
%                also be specified tspan = [t0 t1 ... tN];
%    -  u        is a structure array of nb_traj inputs. u(i) must be a
%                struct with fields
%                    - params_idx : indicates which parameters are made
%                        time dependant
%                    - tin : times when the input changes
%                    - values : values of the parameters
%
%   Outputs:
%
%    -  Pf       Parameter set augmented with the field traj containing
%                 computed trajectories if the input is a param set
%
% or - trajs     array of trajectories
%


% checks if we have a parameter set or a trajectory

output_trajs = 0;
if ~isstruct(P0)
    
    if (size(P0,1) ~= Sys.DimP)
        if (size(P0,2) == Sys.DimP) % be smart, try transpose in case it works
            P0 = P0';
        else
            error('ComputTraj:S0DimensionError',...
                        'Second argument must be a parameter set or be of dimension Sys.DimP x ?')
        end
    end
    output_trajs = 1;
    
    pts = P0;
    P0 = CreateParamSet(Sys,1);
    P0.pts = pts;
    P0.epsi = ones(1, size(pts,2));
end

% checks for an initialization function
if isfield(Sys, 'init_fun')
    P0 = Sys.init_fun(P0);
end

if isfield(P0, 'init_fun')
    P0 = P0.init_fun(P0);
end

if ~isfield(Sys, 'type')
    Sys.type = '';
end

if strcmp(Sys.type,'traces') % No model
    % If system type is only traces, check consistency of params and pts
    Pf=P0;
    for i = 1:numel(Pf.traj)
        Pf.traj(i).param = Pf.pts(1:Pf.DimP,i)';
    end
elseif isfield(P0, 'traj_to_compute') && ~isempty(P0.traj_to_compute) && ~all(P0.traj_to_compute==1:size(P0.pts,2))
    P0 = SPurge(P0);
    S = Sselect(P0, P0.traj_to_compute);
    S = ComputeTraj(Sys, S, tspan);
    Pf = P0;
    Pf.traj = S.traj;
    Pf.Xf = S.Xf;
    return;
end

switch Sys.type

    case 'Extern'
        model = Sys.name;
        Pf = P0;
        ipts = 1:size(P0.pts,2);
        if (numel(ipts)>1)
            fprintf(['Computing ' num2str(numel(ipts)) ' trajectories of model ' model '\n'...
                     '[             25%%           50%%            75%%               ]\n ']);
            iprog =0;
        end
        
        for i= ipts
            if (numel(ipts)>1)
                while (floor(60*i/numel(ipts))>iprog)
                    fprintf('^');
                    iprog = iprog+1;
                end
            end
            
            if isfield(Sys,'init_u')
                U = Sys.init_u(Sys.ParamList(Sys.DimX-Sys.DimU+1:Sys.DimX), P0.pts(1:Sys.DimP,i), tspan);
                assignin('base','t__',U.t);
                assignin('base', 'u__',U.u);
            end
            
            [traj.time, traj.X] = Sys.sim(Sys, tspan, P0.pts(:,i));
            traj.param = P0.pts(1:P0.DimP,i)';
            Pf.traj(i) = traj;
            Pf.Xf(:,i) = traj.X(:,end);
        end
        
        if (numel(ipts)>1)
            fprintf('\n');
        end

    
    case 'Simulink'
        model = Sys.mdl;
        Pf = P0;
        ipts = 1:size(P0.pts,2);
        if (numel(ipts)>1)
            fprintf(['Computing ' num2str(numel(ipts)) ' trajectories of model ' model '\n'...
                     '[             25%%           50%%            75%%               ]\n ']);
            iprog =0;
        end
        
        for i= ipts
            if (numel(ipts)>1)
                while (floor(60*i/numel(ipts))>iprog)
                    fprintf('^');
                    iprog = iprog+1;
                end
            end
            
            if isfield(Sys,'init_u')
                U = Sys.init_u(Sys.ParamList(Sys.DimX-Sys.DimU+1:Sys.DimX), P0.pts(1:Sys.DimP,i), tspan);
                assignin('base','t__',U.t);
                assignin('base', 'u__',U.u);
            end
            
            [traj.time, traj.X] = Sys.sim(Sys, tspan, P0.pts(:,i));
            traj.param = P0.pts(1:P0.DimP,i)';
            Pf.traj(i) = traj;
            Pf.Xf(:,i) = traj.X(:,end);
        end
        
        if (numel(ipts)>1)
            fprintf('\n');
        end
        if ~isfield(Pf, 'traj_ref') % create field traj_ref (one to one mapping)
            Pf.traj_ref = 1:numel(Pf.traj);
        end  
        
    otherwise
        
        InitSystem(Sys);
        
        if iscell(tspan)
            if (numel(tspan)==2)
                T = [tspan{1} tspan{2} tspan{2}];  % not really nice.. should be
                % changed some day
            else
                T= cell2mat(tspan);
            end
        else
            T = tspan;
        end
        
        if (exist('u','var'))
            err = check_u(u);
            if (err~=0)
                error('ComputTraj:ErrorWithU',err);
            end
            
            %This is quite ugly...
            Pf = P0;
            Pf.pts = P0.pts(1:P0.DimP, :);
            Pf=cvm(61, Pf, T, u);
            Pf.pts = P0.pts;
            
            Pf.u = u;
            
        else
            Pf = P0;
            Pf.pts = P0.pts(1:P0.DimP, :);
            Pf=cvm(61, Pf, T); % <- NM: I would love to know how it works inside!
            Pf.pts = P0.pts;
            
  
        end
        
        CVodeFree();
        
        if output_trajs
            Pf = Pf.traj;
        end
end

if isfield(Sys, 'time_mult')
    Pf.time_mult = Sys.time_mult;
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
