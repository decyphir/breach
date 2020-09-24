function [x, fval, output, options] = gbnm(fun,xmin,xmax,options)
%gbnm Globalized and bounded Nelder-Mead as published by Luersen and Le Riche
% 
% Minimizes the scalar function fun in the region determined by
% the column vectors xmin and xmax using a globalized Nelder-Mead
% direct search method based on geometric operations on an N-dimensional
% simplex. 
%
% Usage
%    [x, fval, output, options] = gbnm(fun,xmin,xmax[, options])
%    
% Parameters
%    fun        A handle of the function to be minimized
%    xmin       Column vector of lower bounds
%    xmax       Column vector of upper bounds
%    options    (optional) Structure with optimization algorithm parameters
%
% Returns
%    x          The found global minimum
%    fval       Function value at the found best minimum
%    output     Structure with all found local minima and initial simplices
%    options    Structure of parameters which were used during optimization
%
% Description
%   The method is (almost) identical to the one described in
% M. A. Luersen and R. Le Riche. Globalized nelder-mead method for engineering
% optimization. In ICECT'03: Proceedings of the third international conference on
% Engineering computational technology, pages 165-166, Edinburgh, UK, UK, 2004.
% Civil-Comp press. ISBN 0-948749-84-9
%
% The options structure has the following members and their default values.
% For further explanations of these constants, refer to the source code or to 
% the research paper http://www.emse.fr/~leriche/gbnm_cas.pdf
%       options.maxRestarts = 15; % maximum (probablistic or degenerated) restarts
%       options.maxEvals = 2500;  % maximum function evaluations
%       options.nPoints = 5;      % number of random points per restart
%
%       options.maxIter=250;      % maximum iterations per restart
%       options.alpha = 1;        % reflection coeff
%       options.beta = 0.5;       % contraction coeff
%       options.gamma = 2;        % expansion coeff
%       options.epsilon = 1e-9;   % T2 convergence test coefficient
%       options.ssigma = 5e-4;    % small simplex convergence test coefficient
%     
%
% What else
%   Version    v0.1.1 [2015-09-24] GitHub release
%   Author     Johannes Dorfner
%   License    GPLv3




%--- INIT ---

    if(nargin > 4 || nargin == 0)
        error('Wrong number of arguments for gbnm(fun,xmin,xmax[, options])');
    end

    if nargin < 4
        % Globalized parameters
        options.maxRestarts = 15; % maximum (probablistic or degenerated) restarts
        options.maxEvals = 2500; % maximum function evaluations
        options.nPoints = 5; % number of random points per restart

        % Nelder-mead parameters
        options.maxIter=250; % maximum iterations per restart
        options.alpha = 1; % reflection coeff
        options.beta = 0.5; % contraction coeff
        options.gamma = 2; % expansion coeff
        options.epsilon = 1e-9; % T2 convergence test coefficient
        options.ssigma = 5e-4; % small simplex convergence test coefficient
    end
    
%--- MAIN ---
    
% Pre-allocation
    nEval = 0; 
    Ndim = length(xmin); % number of dimensions
    xrange = xmax - xmin; % range vector
    %p = ones(nProb*ones(1,Ndim))/nProb^Ndim; % probability distribution for restart
    usedPoints = zeros(Ndim, 2*options.maxRestarts);
    usedVals   = Inf*ones(1, 2*options.maxRestarts);
    usedSimplex = cell(1, 2*options.maxRestarts);
    reason = cell(1, options.maxRestarts);
    
    for iRestart=1:options.maxRestarts    
        if (iRestart==1) 
            initial_point = options.xinit;
        else
            % probabilistic restart
            initial_point = probabilistic_restart(usedPoints,xmin,xrange,options.nPoints);
        end
        % probabilistic restart
        %initial_point = probabilistic_restart(usedPoints,xmin,xrange,options.nPoints);
        a=(0.02+0.08*rand)*min(xrange); % simplex size between 2%-10% of min(xrange)
        [Splx, fVal, nEval] = init_simplex(initial_point, a, fun, nEval, xmin, xmax);

        % save initial_point and initial simplex
        usedPoints(:,2*iRestart-1) = Splx(:,1);
        usedVals(2*iRestart-1)     = fVal(1);
        usedSimplex{2*iRestart-1}  = [Splx;fVal];
        
        % simplex iteration
        reason{iRestart} = 'hit maxIter'; % set exitflag
        for iIter=1:options.maxIter
            % sort simplex so that best point is in first column and worst point is last
            [fVal, sortorder] = sort(fVal);
            Splx = Splx(:,sortorder);

            x_worst = Splx(:,end);
            x_best  = Splx(:,1);
            x_cent  = sum(Splx(:,1:Ndim),2)/Ndim; % centroid of all but the worst point

            % T2 converge test
            if std(fVal) < options.epsilon
                reason{iRestart}='convergence (fVals similar)'; 
                break; 
            end
            % small test
            if max((max(Splx,[],2)-min(Splx,[],2))./xrange) < options.ssigma
                reason{iRestart}='convergence (simplex small)';
                break;
            end
            
            if nEval >= options.maxEvals, reason{iRestart}='hit maxEvals'; break; end
            
            % simplex iteration         
            xr = x_cent+options.alpha*(x_cent-x_worst);
            xr = max(xmin,min(xr,xmax)); % respect bounds       
            fr = fun(xr); nEval=nEval+1;
            % see figure 6 in [Luersen2004] for flow diagram
            if fr < fVal(1)
                % reflection better than best point, so try expansion
                xe = x_cent+options.gamma*(xr-x_cent);
                xe = max(xmin,min(xe,xmax)); % respect bounds
                fe = fun(xe); nEval=nEval+1;

                if fe < fr
                    % expansion better than reflection, so use xe
                    Splx(:,end) = xe;
                    fVal(end) = fe;
                else
                    % expansion worse than reflection, so use xr
                    Splx(:,end) = xr;
                    fVal(end) = fr;
                end
                
            else
                % reflection not better than best point
                if fr <= fVal(end-1)
                    % reflection still better than second worst point, so use xr
                    Splx(:,end) = xr;
                    fVal(end) = fr;
                else
                    % reflection worse than second worst point, so contract
                    if fr < fVal(end)
                        % but still use xr if it is better than worst point
                        Splx(:,end) = xr;
                        fVal(end) = fr;
                    end
                    % contraction of worst point towards centroid
                    xc = x_cent+options.beta*(x_worst-x_cent);
                    fc = fun(xc); nEval=nEval+1;

                    if fc > fVal(end)
                        % xc still worse than worst point, so contract all
                        % simplex points towards best point (NOT centroid!)
                        for k=2:Ndim+1
                            Splx(:,k)=0.5*(Splx(:,k) + x_best);
                            fVal(k) = fun(Splx(:,k)); nEval=nEval+1;
                        end
                    else
                        % xc at least as good as worst point, so use it
                        Splx(:,end) = xc;
                        fVal(end) = fc;
                    end
                end
            end
        end % for iIter
        
        % store best point
        usedPoints(:,2*iRestart) = Splx(:,1);
        usedVals(2*iRestart)     = fVal(1);
        usedSimplex{2*iRestart}  = [Splx;fVal];
        
        % output progress
%         disp([num2str(iRestart) ':']);
%         disp(['   Reason:  ' reason]);
%         disp(['   Vertex:  ' num2str(Splx(:,1)')]);
%         disp(['   Value:   ' num2str(fVal(1))]);
%         disp(['   nEval:   ' num2str(nEval) ' (' num2str(options.maxEvals) ')']);
%         drawnow;
        
        % if maxEvals is hurt also break from restart loop
        if nEval >= options.maxEvals, reason{iRestart}='hit maxEvals'; break; end   
    end % for iRestart

    
    
    % finally output globally best solution (and all other found optima)
    [fval, minidx] = min(usedVals);
    x = usedPoints(:,minidx);
    % output convergence points, then starting points
    usedPoints = usedPoints(:,[2:2:end 1:2:end]);
    usedVals = usedVals(:,[2:2:end 1:2:end]);

    output = struct(...
        'usedPoints',usedPoints,...
        'usedVals',usedVals,...
        'usedSimplex',{usedSimplex},...
        'reason',{reason},...
        'nEval',nEval...
        );
end % function gbnm

% -------------------------------------------------------------------------




% SUBFUNCTIONS

function [Splx, fVal, nEval] = init_simplex(x0, a, fun, nEval, xmin, xmax)
%init_simplex initialises new simplex at point x0 with size a
% fun is used to calculate fVal at all points
    Ndim=length(x0);
    Splx=zeros(Ndim,Ndim+1);
    fVal=zeros(1,Ndim+1);
    
    p=a*(sqrt(Ndim+1)+Ndim-1)/(Ndim*sqrt(2));
    q=a*(sqrt(Ndim+1)     -1)/(Ndim*sqrt(2));
    
    Splx(:,1)=x0;
    for k=1:Ndim 
        Splx(:,k+1)=Splx(:,1);
        Splx(k,k+1)=Splx(k,k+1)+p;
        Splx([1:k-1 k+1:Ndim],k+1)=Splx([1:k-1 k+1:Ndim],k+1)+q*ones(Ndim-1,1);
        
        Splx(:,k+1)=max(xmin,min(Splx(:,k+1),xmax)); % respect bounds
    end
    % evaluate function for all points
    for k=1:Ndim+1
        fVal(k)= fun(Splx(:,k));
        nEval=nEval+1;
    end
end


function best_point = probabilistic_restart(usedPoints,xmin,xrange,nPoints)
    Ndim=size(usedPoints,1);
    best_prob=Inf;
    for k=1:nPoints
        % cast a random point somewhere in the allowed range     
        random_point=xmin+xrange.*rand(Ndim,1);
        random_prob = gauss(random_point, usedPoints, xrange);

        if random_prob < best_prob % minimise probability
            best_point = random_point;
            best_prob  = random_prob;
        end
    end
end


function prob = gauss(x,points,xrange)
    glp = 0.01; % Gaussian length parameter (see [Luersen2004, 2.1 Probabilistic restart])
    [Ndim, Npts] = size(points);
    sigma = diag(glp*xrange.^2);
    sigmainv = diag(1./diag(sigma));
    prob = 0;
    for k=1:Npts
        mu = points(:,k);
        if mu == zeros(Ndim,1), break; end
        prob = prob + exp(-.5*(x-mu)'*sigmainv*(x-mu))/sqrt((2*pi)^Ndim*det(sigma));
    end
end