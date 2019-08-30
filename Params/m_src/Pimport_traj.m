function P = Pimport_traj(P, P0)
%Pimport_traj get trajectories from P0 in P, update traj_ref and traj_to_compute 
%

    if isfield(P0, 'traj')
        if isfield(P,'traj')
            P.traj = [P.traj P0.traj];
        else
            P.traj = P0.traj;          
        end
        P = Preset_traj_ref(P);
    end
    
