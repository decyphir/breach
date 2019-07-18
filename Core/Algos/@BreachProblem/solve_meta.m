function res = solve_meta(this)

% display header
if ~strcmp(this.display,'off')
    fprintf('\nSTART OPTIMIZATION METAHEURISTICS\n');   
end

opt = this.solver_options;

if ~strcmp(this.display,'off')
    fprintf('\n     TEST CORNERS\n');
end

if opt.num_corners>0
    res= this.solve_corners();
end

if ~strcmp(this.display,'off')&&~this.stopping()
    this.display_status();    
    this.Display_Best_Results(res.f, res.x);
end

while ~this.stopping()
    qr_step = opt.quasi_rand_seed;
    
    %% Quasi-random phase
    if ~strcmp(this.display,'off')
        fprintf('\nTEST QUASI-RANDOM SAMPLES\n');
    end

    this.setup_quasi_random('quasi_rand_seed',qr_step,'num_quasi_rand_samples',opt.num_quasi_rand_samples);
    res = this.solve_quasi_random();
    opt.quasi_rand_seed = opt.quasi_rand_seed+opt.num_quasi_rand_samples;  

    %% Disp
    if this.stopping()
        break;
    else
        if ~strcmp(this.display,'off')
            this.Display_Best_Results(res.fval, res.x);
        end
    end
    x0 = res.x; % best x for next phase
    
    %% Local phase
    if ~strcmp(this.display,'off')
        fprintf('\nRUN LOCAL OPTIMIZATION\n');
    end

    res = run_nelder_mead(this,opt, x0); 
    this.solver_options = opt;    
    %% Disp 
    if this.stopping()        
        break;
    else
        if ~strcmp(this.display,'off')
            this.Display_Best_Results(res.fval, res.x);
        end
    end
end

if ~strcmp(this.display,'off')
    fprintf('\nEND OPTIMIZATION METAHEURISTICS\n');
end

end

function res = run_nelder_mead(this, opt,  x0)

this.setup_nelder_mead(x0, 'MaxFunEvals', opt.local_max_obj_eval, 'Display', 'off');
res = this.solve_nelder_mead();    

end



