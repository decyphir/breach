function U = SimulinkInput(input_opt, pts, tspan)

cp = input_opt.cp;
pts_u = pts(input_opt.idx);

switch input_opt.type
    
    case 'UniStep'
        U.u = [];
        U.t = tspan';
        for i_cp = 1:numel(input_opt.cp)
            cp_values = pts_u(1:cp(i_cp));
            method = input_opt.method{i_cp};
            pts_u = pts_u(cp(i_cp)+1:end);
            Ui = GetUi(method);
            U.u = [U.u, Ui.u];
        end
        
    case 'VarStep'
        U.u = [];
        U.t = tspan';
        for i_cp = 1:numel(input_opt.cp)
            method = input_opt.method{i_cp};
            cp_values = pts_u(1:2*cp(i_cp));
            pts_u = pts_u(2*cp(i_cp)+1:end);
            Ui = GetVarUi(method);
            U.u = [U.u, Ui.u];
        end
        
end

    function Ui = GetUi(method)
        
        Ui.t = tspan';
        t_cp = linspace(tspan(1), tspan(end), cp(i_cp)+1)';
        if numel(t_cp)==2
            Ui.u = cp_values(1)*ones(numel(tspan),1);
        else
            Ui.u = interp1(t_cp(1:end-1), cp_values, tspan', method, 'extrap');
        end
    end

    function Ui = GetVarUi(method)
        
        Ui.t = tspan';
        dt_cp = cp_values(1:2:end);
        t_cp = unique( [0; cumsum(dt_cp)]);
        u_values = cp_values(2:2:end);
        u_values = u_values(1:min([numel(t_cp) numel(u_values)]));
        t_cp = t_cp(1:min([numel(t_cp) numel(u_values)]));
        if numel(t_cp)==1
            Ui.u = u_values(1)*ones(numel(tspan),1);
        else
            Ui.u = interp1(t_cp, u_values, tspan', method, 'extrap');
        end
        
    end


end