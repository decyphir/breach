function res = STL_EvalTemplate(Sys, phi, P, trajs, sigs_phi, sigs_sys)
%STL_EVALTEMPLATE computes the satisfaction function of a property for one or
% many traces.
%
% Synopsis: res = STL_EvalTemplate(Sys, phi, phi_sigs, P, trajs[, sys_sigs])
%
% Inputs:
%  - Sys        : the system
%  - phi        : a STL Formula
%  - P          : is a parameter set which contains one parameter vector only
%                 used for properties parameters.
%  - trajs      : is a structure with fields X and time. It may contains many
%                 trajectories. In this case, all will be checked wrt the
%                 property parameter described in P.
%  - sigs_phi   : templated signals
%  - sigs_sys   : list of signals to check (default: all)
%
% Outputs: a structure res with fields
%  - all  : list of signals satisfying the templates for all trajs
%  - some : list of pairs of signal and traj number satisfying template
%  - none : list of signals satisfying the templates for no trajs
%

%% Arguments

%% Switch on number of template signals, start with 1 and 2 first

ParamList = Sys.ParamList;

res.all  = {};
res.some = {};
res.none = {};

id_phi = get_id(phi);

switch numel(sigs_phi)
    case 1
        
        if ~exist('sigs_sys','var')
            sigs_sys = 1:Sys.DimX;
        end
        
        if iscell(sigs_sys)
            idx_sig_sys=  FindParam(Sys,sigs_sys);
        else
            idx_sig_sys = sigs_sys;
        end
        
        for i_sig_sys = idx_sig_sys
            sig_sys=  ParamList{i_sig_sys};
            Sys.ParamList = ParamList;
            P.ParamList = ParamList;
            Sys.ParamList{i_sig_sys} = sigs_phi{1};
            P.ParamList{i_sig_sys} = sigs_phi{1};
            val =  STL_Eval(Sys, phi, P, trajs, 0);
%            fprintf('Checking %s for %s=%s', id_phi,sigs_phi{1}, sig_sys);
%            fprintf('---> rob=%g \n', val);

            if all(val>0)
                res.all = [res.all {sig_sys}];
            elseif any(val>0)
                res.some = [res.some {sig_sys}];
            else
                res.none = [res.none {sig_sys}];
            end
        end
        
    case 2 % assumes sys_sigs is a cell of cell
        
        if ~exist('sigs_sys','var')
            sigs_sys = { 1:Sys.DimX, 1:Sys.DimX};
        end
        
        if iscell(sigs_sys{1})
            idx_sig_sys1=  FindParam(Sys,sigs_sys{1});
            idx_sig_sys2=  FindParam(Sys,sigs_sys{2});
        else
            idx_sig_sys1 = sigs_sys{1};
            idx_sig_sys2 = sigs_sys{2};
        end
        
        for i_sig_sys1 = idx_sig_sys1
            for i_sig_sys2 = idx_sig_sys2
                if i_sig_sys2~= i_sig_sys1     %  otherwise why have two signals template?
    
                    sig_sys1=  ParamList{i_sig_sys1};
                    sig_sys2=  ParamList{i_sig_sys2};
                    Sys.ParamList = ParamList;
                    P.ParamList = ParamList;
                    
                    Sys.ParamList{i_sig_sys1} = sigs_phi{1};
                    P.ParamList{i_sig_sys1} = sigs_phi{1};
                    
                    Sys.ParamList{i_sig_sys2} = sigs_phi{2};
                    P.ParamList{i_sig_sys2} = sigs_phi{2};
                    val =  STL_Eval(Sys, phi, P, trajs, 0);
                    
                    %fprintf('Checking %s for %s=%s and %s=%s', id_phi,sigs_phi{1}, sig_sys1, sigs_phi{2}, sig_sys2); 
                    %fprintf('---> rob=%g \n', val);

                    if all(val>0)
                        res.all = [res.all {{sig_sys1, sig_sys2}}];
                    elseif any(val>0)
                        res.some = [res.some {{sig_sys1, sig_sys2}}];
                    else
                        res.none = [res.none {{sig_sys1, sig_sys2}}];
                    end
                end
            end
        end
        
        
end
