%MINIMIZE        Solve constrained optimization problems, 
%                globally or locally
%
% Usage:
%    sol = MINIMIZE(func, x0) 
%    sol = MINIMIZE(..., x0, A, b) 
%    sol = MINIMIZE(..., b, Aeq, beq) 
%    sol = MINIMIZE(..., beq, lb, ub)
%    sol = MINIMIZE(..., ub, nonlcon) 
%    sol = MINIMIZE(..., nonlcon, options) 
%
% [sol, fval] = MINIMIZE(func, ...)
% [sol, fval, exitflag] = MINIMIZE(func, ...)
% [sol, fval, exitflag, output] = MINIMIZE(func, ...)
% [sol, fval, exitflag, output, grad] = MINIMIZE(func, ...)
% [sol, fval, exitflag, output, gradient, hessian] = MINIMIZE(func, ...)
%
% INPUT ARGUMENTS:
%
%  fun, x0 - see FMINSEARCH or FMINLBFGS.
%
%  A, b - (OPTIONAL) Linear inequality constraint array and right
%       hand side vector. 
%
%  Aeq, beq - (OPTIONAL) Linear equality constraint array and right
%       hand side vector. 
%
%  lb, ub - (OPTIONAL) lower/upper bound vector or array. Both must have 
%       the same size as x0.
%
%       If no lower bounds exist for one of the variables, then
%       supply -inf for that variable. Similarly, if no upper bounds 
%       exist, supply +inf. If no bounds exist at all, then [lb] and/or 
%       [ub] may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
%  nonlcon - (OPTIONAL) function handle to general nonlinear constraints,
%       inequality and/or equality constraints.
%
%       [nonlcon] must return two vectors, [c] and [ceq], containing the
%       values for the nonlinear inequality constraints [c] and
%       those for the nonlinear equality constraints [ceq] at [x]. (Note: 
%       these constraints were chosen to be consistent with those of 
%       fmincon.)
%
%  options - (OPTIONAL) an options structure created manually or with
%       setoptimoptions().    
%
% OUTPUT ARGUMENTS:
%
% sol, fval - the solution vector and the corresponding function value,
%       respectively. 
%
% exitflag - (See also the help on FMINSEARCH) A flag that specifies the 
%       reason the algorithm terminated. FMINSEARCH uses only the values
%
%           1    fminsearch converged to a solution x
%           0    Max. # of function evaluations or iterations exceeded
%          -1    Algorithm was terminated by the output function.
%
%       Since MINIMIZE handles constrained problems, the following 
%       values were added:
%
%           2    Problem overconstrained by either [lb]/[ub] or
%                [Aeq]/[beq] - nothing done
%          -2    Problem is infeasible after the optimization (some or 
%                any of the constraints are violated at the final 
%                solution).
%          -3    INF or NAN encountered during the optimization. 
%
% output - (See also the help on FMINSEARCH) A structure that contains
%       additional details on the optimization. 
%
% gradient/hessian - (only for FMINLBFGS) the gradient and estimated 
%                    hessian at x = sol.
%
% Notes:
%
%  If [options] is supplied, then TolX will apply to the transformed
%  variables. All other FMINSEARCH parameters should be unaffected.
%
%
% EXAMPLES:
%
% rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%
% >> % Fully unconstrained problem
% >> minimize(rosen, [3 3])
% ans =
%    1.0000    1.0000
%
%
% >> % lower bound constrained
% >> minimize(rosen,[3 3], [],[], [],[], [2 2])
% ans =
%    2.0000    4.0000
%
%
% >> % x(2) fixed at 3
% >> minimize(rosen,[3 3], [],[], [],[], [-inf 3],[inf,3])
% ans =
%    1.7314    3.0000
%
%
% >> % simple linear inequality: x(1) + x(2) <= 1
% >> minimize(rosen,[0; 0], [1 1], 1)
% 
% ans =
%    0.6187    0.3813
% 
%
% >> % nonlinear inequality: sqrt(x(1)^2 + x(2)^2) <= 1
% >> % nonlinear equality  : x(1)^2 + x(2)^3 = 0.5
%
% Execute this m-file:
%
%   function test_minimize
%        rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%
%        options = optimset('TolFun', 1e-8, 'TolX', 1e-8);
%
%        minimize(rosen, [3 3], [],[],[],[],[],[],...
%        @nonlcon, [], options)
%
%   end
%   function [c, ceq] = nonlcon(x)
%        c = norm(x) - 1;
%        ceq = x(1)^2 + x(2)^3 - 0.5;
%   end
%
% ans =
%    0.6513    0.4233
%
%
% Of course, any combination of the above constraints is

% also possible.
%
%
% See also: SETOPTIMOPTIONS, FMINSEARCH, FMINLBFGS.


% Please report bugs and inquiries to: 
%
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com         (personal)
%              r-oldenhuis@ispace-inc.com  (professional)
% Affiliation: ispace inc., Japan
% License    : BSD 2.0; see license.txt


% FMINSEARCHBND, FMINSEARCHCON and part of the documentation 
% for MINIMIZE witten by
%
% Author : John D'Errico
% E-mail : woodchips@rochester.rr.com


% TODO
%{
WISH: Transformations similar to the ones used for bound constraints should 
      also be possible for linear constraints; figure this out.
WISH: more checks for inconsistent constraints: 
      - equality constraints might all lie outside the region defined by the 
        bounds and linear inequalities
      - different linear equality constraints might never intersect
WISH: Include examples and demos in the documentation
WISH: a properly working ConstraintsInObjectiveFunction

FIXME: hess not ouput AT ALL
FIXME: check how functions are evaluated; isn't it better to have objFcn()
       and conFcn() do more work? 
FIXME: ignore given nonlcon() if ConstraintsInObjectiveFunction is true?
FIXME: fevals is miscounted; possibly fminlbfgs bug
FIXME: TolX, TolFun, TolCon, DiffMinchange, DiffMaxChange should be transformed
FIXME: 'AlwaysHonorConstraints' == 'bounds' does nothing
%}
function [sol, fval,...
          exitflag, output,...
          grad, hess] = minimize(funfcn, x0,...
                                 A,b,...
                                 Aeq,beq,...
                                 lb,ub,...
                                 nonlcon,...
                                 options, varargin)
             
    % If you find this work useful, please consider a donation:
    % https://www.paypal.me/RodyO/3.5
        
    %% Initialization
    
    % Process user input
    narg  = nargin;
    nargo = nargout;
    if verLessThan('MATLAB', '8.6')
        error(nargchk(2, inf, narg,    'struct'));  %#ok<NCHKN>
        error(nargchk(0,   6, nargout, 'struct'));  %#ok<NCHKM>
    else
        narginchk (2, inf);
        nargoutchk(0, 6);
    end
    
    if (narg < 10) || isempty(options), options = setoptimoptions; end
    if (narg <  9), nonlcon = ''; end
    if (narg <  8), ub  = []; end
    if (narg <  7), lb  = []; end
    if (narg <  6), beq = []; end
    if (narg <  5), Aeq = []; end
    if (narg <  4), b   = []; end
    if (narg <  3), A   = []; end
            
    % Extract options
    nonlconFcn_in_objFcn = getoptimoptions('ConstraintsInObjectiveFunction', false);
    tolCon        = getoptimoptions('TolCon', 1e-8);   
    algorithm     = getoptimoptions('Algorithm', 'fminsearch');
    strictness    = getoptimoptions('AlwaysHonorConstraints', 'none');    
    diffMinChange = getoptimoptions('DiffMinChange', 1e-8);
    diffMaxChange = getoptimoptions('DiffMaxChange', 1e-1);        
    finDiffType   = getoptimoptions('FinDiffType', 'forward');        
    OutputFcn     = getoptimoptions('OutputFcn', []);
    PlotFcn       = getoptimoptions('PlotFcn', []);
            
    % Set some logicals for easier reading
    have_nonlconFcn    = ~isempty(nonlcon) || nonlconFcn_in_objFcn;
    have_linineqconFcn = ~isempty(A) && ~isempty(b);
    have_lineqconFcn   = ~isempty(Aeq) && ~isempty(beq); 
    
    create_output = (nargout >= 4);
    
    need_grad                    = strcmpi(algorithm,'fminlbfgs');
    grad_obj_from_objFcn         = need_grad && strcmpi(options.GradObj   ,'on');
    grad_nonlcon_from_nonlconFcn = need_grad && strcmpi(options.GradConstr,'on') && have_nonlconFcn;
    
    do_display          = ~isempty(options.Display) && ~strcmpi(options.Display,'off');
    do_extended_display = do_display && strcmpi(options.Display,'iter-detailed');
    
    do_global_opt = isempty(x0);
    
    % Do we have an output function? 
    have_outputFcn = ~isempty(OutputFcn);
    have_plotFcn   = ~isempty(PlotFcn);
        
    % No x0 given means minimize globally.
    if do_global_opt
        [sol, fval, exitflag, output, grad, hess] = minimize_globally();
        return;
    end
            
    % Make copy of UNpenalized function value
    UPfval = inf;    
    
    % Initialize & define constants    
    superstrict = false;      % initially, don't use superstrict setting    
    exp50   = exp(50);        % maximum penalty
    N0      = numel(x0);      % variable to check sizes etc.
    Nzero   = zeros(N0, 1);   % often-used zero-matrix    
    grad    = [];             % initially, nothing for gradient
    sumAT   = repmat(sum(A,1).',1,size(b,2));     % column-sum of [A]  , transposed and replicated
    sumAeqT = repmat(sum(Aeq,1).',1,size(beq,2)); % column-sum of [Aeq], transposed and replicated
    
    % Initialize output structure
    output = [];
    if create_output
        
        % Fields always present
        output.iterations = 0;
        output.algorithm  = '';
        output.message    = 'Initializing optimization...';
        
        % Fields present depending on the presence of nonlcon
        if ~have_nonlconFcn
            output.funcCount = 0;
        else
            output.ObjfuncCount = 0;
            output.ConstrfuncCount = 1; % one evaluation in check_input()
            output.constrviolation.nonlin_eq   = cell(2,1);
            output.constrviolation.nonlin_ineq = cell(2,1);
        end
        
        % Fields present depending on the presence of A or Aeq
        if have_linineqconFcn
            output.constrviolation.lin_ineq = cell(2,1); end
        if have_lineqconFcn
            output.constrviolation.lin_eq = cell(2,1); end
    end
    
    % Save variable with original size 
    new_x = x0;
        
    % Check for an output/plot functions. If there are any, 
    % use wrapper functions to call it with un-transformed variable
    if have_outputFcn
        OutputFcn = options.OutputFcn;
        if ~iscell(OutputFcn)
            OutputFcn = {OutputFcn}; end
        options.OutputFcn = @OutputFcn_wrapper;
    end    
    if have_plotFcn
        PlotFcn = options.PlotFcn;
        if ~iscell(PlotFcn)
            PlotFcn = {PlotFcn}; end
        options.PlotFcn = @PlotFcn_wrapper;
    end  
    
    % Adjust bounds when they are empty
    if isempty(lb), lb = -inf(size(x0)); end
    if isempty(ub), ub = +inf(size(x0)); end
           
    % Check the user-provided input with nested function check_input    
    if ~check_input()
        return; end        
        
    % Force everything to be column vector
    ub = ub(:);   x0    = x0(:);    
    lb = lb(:);   x0one = ones(size(x0));           
    
    % replicate lb or ub when they are scalars, and x0 is not
    if isscalar(lb) && (N0 ~= 1), lb = lb*x0one; end
    if isscalar(ub) && (N0 ~= 1), ub = ub*x0one; end    
    
    % Determine the type of bounds
    nf_lb   = ~isfinite(lb);                nf_ub   = ~isfinite(ub); 
    fix_var =  lb == ub;                    lb_only = ~nf_lb &  nf_ub & ~fix_var;
    ub_only =  nf_lb & ~nf_ub & ~fix_var;   unconst =  nf_lb &  nf_ub & ~fix_var;
    lb_ub   = ~nf_lb & ~nf_ub & ~fix_var;    
    
    
    %% Optimization 
    
    % Force the initial estimate inside the given bounds
    x0(x0 < lb) = lb(x0 < lb);  x0(x0 > ub) = ub(x0 > ub);
        
    % Transform initial estimate to its unconstrained counterpart
    xin = x0;                                        % fixed and unconstrained variables   
    xin(lb_only) = sqrt(x0(lb_only) - lb(lb_only));  % lower bounds only   
    xin(ub_only) = sqrt(ub(ub_only) - x0(ub_only));  % upper bounds only
    xin(lb_ub)   = real(asin( 2*(x0(lb_ub) - lb(lb_ub))./ ...
                     (ub(lb_ub) - lb(lb_ub)) - 1));  % both upper and lower bounds
    xin(fix_var) = [];
    
    % Some more often-used matrices
    None    = ones(numel(xin)+1,1);  
    Np1zero = zeros(N0, numel(xin)+1);
   
    % Optimize the problem
    try
        switch lower(algorithm)

            % MATLAB's own derivative-free Nelder-Mead algorithm (FMINSEARCH())
            case 'fminsearch'
                [presol, fval, exitflag, output_a] = ...
                    fminsearch(@funfcnT, xin, options);                          

                % Transform solution back to original (bounded) variables...
                sol = new_x;    sol(:) = X(presol); % with the same size as the original x0

                % Evaluate function once more to get unconstrained values
                % NOTE: this eval is added to total fevals later
                fval(:) = objFcn(sol);

            % Steepest descent or Quasi-Newton (limited-memory) BFGS 
            % (both using gradient) FMINLBFGS(), by Dirk-Jan Kroon 
            case 'fminlbfgs'

                % DEBUG: check gradients
                %{

                % Numerical
                oneHundredth = [
                    ( funfcnT(xin+[1e-2;0])-funfcnT(xin-[1e-2;0]))/2e-2
                    ( funfcnT(xin+[0;1e-2])-funfcnT(xin-[0;1e-2]))/2e-2]

                oneTrillionth = [
                    ( funfcnT(xin+[1e-12;0])-funfcnT(xin-[1e-12;0]))/2e-12
                    ( funfcnT(xin+[0;1e-12])-funfcnT(xin-[0;1e-12]))/2e-12]

                h = 1e-3;
                dfdx   = @(f,x) ( (f(x-4*h)-f(x+4*h))/280 + 4*(f(x+3*h)-f(x-3*h))/105 + (f(x-2*h)-f(x+2*h))/5 + 4*(f(x+h)-f(x-h))/5 )/h;  %#ok
                highOrder = [ dfdx(@(x)funfcnT([x;xin(2)]), xin(1)); dfdx(@(x)funfcnT([xin(1);x]), xin(2)) ]

                % As computed in penalized/transformed function
                 [F,G] = funfcnT(xin)
                %}

                [presol, fval, exitflag, output_a, ~, hess_a] = ...
                    fminlbfgs(@funfcnT, xin, options);

                % Transform solution back to original (bounded) variables
                sol    = new_x;   
                sol(:) = X(presol); % with the same size as the original x0  
                
                % Evaluate function some more to get unconstrained values            
                if grad_obj_from_objFcn
                    % Function value and gradient
                    [fval, grad] = objFcn(sol);
                    if create_output
                        output_a.funcCount = output_a.funcCount + 1; end
                else
                    % Only function value
                    [fevals, grad] = computeJacobian(funfcn, sol, objFcn(sol));
                    if create_output
                        output_a.funcCount = output_a.funcCount + fevals + 1; end
                end

        end % switch (algorithm)
        
    catch ME                
        ME2 = MException([mfilename ':unhandled_error'], [...
                         'Unhandled problem encountered; please contact ',...
                         'the author with this exact message, and your ',...
                         'exact inputs.']);               
        throw(addCause(ME2,ME));
    end
    
    % Copy appropriate fields to the output structure
    if create_output
        output.message    = output_a.message;
        output.algorithm  = output_a.algorithm;
        output.iterations = output_a.iterations;
        if ~have_nonlconFcn
            output.funcCount    = output_a.funcCount + 1;
        else
            output.ObjfuncCount = output_a.funcCount + 1;
        end
    end
    
    % Append constraint violations to the output structure, and change the
    % exitflag accordingly
    [output, exitflag] = finalize(sol, output, exitflag);
    
    % Hessian output
    if strcmpi(algorithm, 'fminlbfgs')    
        
        % No violation - Hessian estimate equals last-evaluated
        if exitflag ~= -2
            hess = hess_a;
            
        % Constraints are violated - no unconstrained estimate is available
        else
            hess = NaN(size(hess_a));
        end
    end
    
    %% NESTED FUNCTIONS (THE ACTUAL WORK)
    
    % Check user-provided input
    function go_on = check_input()
        
% FIXME: diffMinChange and diffMaxChange can be inconsistent
        
        go_on = true;
        
        % Dimensions & weird input
        if (numel(lb) ~= N0 && ~isscalar(lb)) || (numel(ub) ~= N0 && ~isscalar(ub))
            error([mfilename ':lb_ub_incompatible_size'], [...
                  'Size of either [lb] or [ub] incompatible with size ',...
                  'of [x0].'])
        end
        if ~isempty(A) && isempty(b)
            warning([mfilename 'minimize:Aeq_but_not_beq'], [...
                    'I received the matrix [A], but you omitted the ',...
                    'corresponding vector [b].\nI''ll assume a zero-vector ',...
                    ' for [b]...']);
            b = zeros(size(A,1), size(x0,2));
        end
        if ~isempty(Aeq) && isempty(beq)
            warning([mfilename ':Aeq_but_not_beq'], [...
                    'I received the matrix [Aeq], but you omitted the ',...
                    'corresponding vector [beq].\nI''ll assume a ',...
                    'zero-vector for [beq]...']);
            beq = zeros(size(Aeq,1), size(x0,2));
        end
        if isempty(Aeq) && ~isempty(beq)
            warning([mfilename ':beq_but_not_Aeq'],[...
                    'I received the vector [beq], but you omitted the ',...
                    'corresponding matrix [Aeq].\nI''ll ignore the given ',...
                    '[beq]...']);
            beq = [];
        end
        if isempty(A) && ~isempty(b)
            warning([mfilename ':b_but_not_A'], [...
                    'I received the vector [b], but you omitted the ',...
                    'corresponding matrix [A].\nI''ll ignore the given ',...
                    '[b]...']);
            b = [];
        end
        if have_linineqconFcn && size(b,1)~=size(A,1)
            error([mfilename ':b_incompatible_with_A'],...
                  'The size of [b] is incompatible with that of [A].');
        end
        if have_lineqconFcn && size(beq,1)~=size(Aeq,1)
            error([mfilename ':b_incompatible_with_A'],...
                  'The size of [beq] is incompatible with that of [Aeq].');
        end
        if ~isvector(x0) && ~isempty(A) && (size(A,2) ~= size(x0,1))
            error([mfilename ':A_incompatible_size'], [...
                  'Linear constraint matrix [A] has incompatible size for ',... 
                  'given [x0].']);
        end
        if ~isvector(x0) && ~isempty(Aeq) && (size(Aeq,2) ~= size(x0,1))
            error([mfilename ':Aeq_incompatible_size'], [...
                  'Linear constraint matrix [Aeq] has incompatible size ',...
                  'for given [x0].']);
        end
        if ~isempty(b) && size(b,2)~=size(x0,2)
            error([mfilename ':x0_vector_but_not_b'],[...
                  'Given linear constraint vector [b] has incompatible ',...
                  'size with given [x0].']);
        end
        if ~isempty(beq) && size(beq,2)~=size(x0,2)
            error([mfilename ':x0_vector_but_not_beq'], [...
                  'Given linear constraint vector [beq] has incompatible ',...
                  'size with given [x0].']);
        end 
              
        % Functions are not function handles        
        if ~isa(funfcn, 'function_handle')
            error([mfilename ':func_not_a_function'],...
                  'Objective function must be given as a function handle.');
        end
        if ~isempty(nonlcon) && ~ischar(nonlcon) && ~isa(nonlcon, 'function_handle')
            error([mfilename ':nonlcon_not_a_function'], [...
                  'Non-linear constraint function must be a function ',...
                  'handle (advised) or string (discouraged).']);
        end
        
        % Check if FMINLBFGS can be executed
        if ~isempty(algorithm) && strcmpi(algorithm,'fminlbfgs') 
            assert(~isempty(which('fminlbfgs')),...
                   [mfilename ':fminlbfgs_not_present'],[...
                   'The function FMINLBFGS is not present in the current ',...
                   'MATLAB path.']);
        end
        
        % check output arguments
        if nargo >= 5
            assert(strcmpi(algorithm,'fminlbfgs'),...
                   [mfilename ':fminlbfgs_needed_for_derivative_output'], [...
                   'Gradient or Hessian are only used by the fminlbfgs ',...
                   'algorithm. Set ''algorithm'' to ''fminlbfgs'' when ',...
                   'requesting derivative information as an output.']);
        end
        

        % evaluate the non-linear constraint function on the initial value, 
        % to perform initial checks
        grad_c   = [];  
        grad_ceq = [];
        [c, ceq] = conFcn(x0);
        incrementConfuncCount();

        % Check sizes of derivatives
        if ~isempty(grad_c) && (size(grad_c,2) ~= numel(x0)) && (size(grad_c,1) ~= numel(x0))
            error([mfilename ':grad_c_incorrect_size'], [...
                  'The matrix of gradients of the non-linear in-equality ',...
                  'constraints\nmust have one of its dimensions equal to ',...
                  'the number of elements in [x].']);            
        end
        if ~isempty(grad_ceq) && (size(grad_ceq,2) ~= numel(x0)) && (size(grad_ceq,1) ~= numel(x0))
             error([mfilename ':grad_ceq_incorrect_size'], [...
                   'The matrix of gradients of the non-linear equality ',...
                   'constraints\nmust have one of its dimension equal to ',...
                   'the number of elements in [x].']);
        end        
        
        % Test the feasibility of the initial solution (when strict or
        % superstrict behavior has been enabled)        
        if strcmpi(strictness, 'Bounds') || strcmpi(strictness, 'All')
            
            superstrict = strcmpi(strictness, 'All');   
            
            if ~isempty(A) && any(any(A*x0 > b))
                error([mfilename ':x0_doesnt_satisfy_linear_ineq'],[...
                      'Initial estimate does not satisfy linear inequality.', ...
                      '\nPlease provide an initial estimate inside the ',...
                      'feasible region.']);
            end
            if ~isempty(Aeq) && any(any(Aeq*x0 ~= beq))
                error([mfilename ':x0_doesnt_satisfy_linear_eq'],[...
                      'Initial estimate does not satisfy linear equality.', ...
                      '\nPlease provide an initial estimate inside the ',...
                      'feasible region.']);
            end
            if have_nonlconFcn
                % check [c]
                if ~isempty(c) && any(c(:) > ~superstrict*tolCon) 
                    error([mfilename ':x0_doesnt_satisfy_nonlinear_ineq'], [...
                          'Initial estimate does not satisfy nonlinear ',...
                          'inequality.\nPlease provide an initial estimate ',...
                          'inside the feasible region.']);
                end
                % check [ceq]
                if ~isempty(ceq) && any(abs(ceq(:)) >= ~superstrict*tolCon) 
                    error([mfilename ':x0_doesnt_satisfy_nonlinear_eq'], [...
                          'Initial estimate does not satisfy nonlinear ',...
                          'equality.\nPlease provide an initial estimate ',...
                          'inside the feasible region.']);
                end
            end
        end    
                
        % Detect and handle degenerate problems
        
% FIXME: with linear constraints it should be easy to determine if the
% constraints make feasible solutions impossible...
        
        % Impossible constraints
        inds = all(A==0,2);
        if any(inds) && any(any(b(inds,:) ~= 0))
            error([mfilename ':impossible_linear_inequality'],...
                  'Impossible linear inequality specified.');
        end
        
        inds = all(Aeq==0,2);
        if any(inds) && any(any(beq(inds,:)~=0))
            error([mfilename ':impossible_linear_equality'],...
                  'Impossible linear equality specified.');
        end
        
        % Degenerate constraints
        if size(Aeq,1) > N0 && rank(Aeq) == N0
            warning([mfilename ':linear_equality_overconstrains'], [...
                    'Linear equalities overconstrain problem; constrain ',...
                    'violation is likely.']);
        end
        
        if size(Aeq,2) >= N0 && rank(Aeq) >= N0
            
            warning([mfilename ':linear_equality_overconstrains'],...
                    'Linear equalities define solution - nothing to do.');
            
            sol   = Aeq\beq;
            fval  = objFcn(sol);
            new_x = sol;
            
            exitflag = 2;
            if create_output
                output.iterations = 0;
                output.message = 'Linear equalities define solution; nothing to do.';
                incrementObjfuncCount();
            end
            
            do_display_P = do_display;
            do_display   = false;            
            [output, exitflag] = finalize(sol, output, exitflag);
            
            if create_output && exitflag ~= -2
                output.message = sprintf(['%s\nFortunately, the solution ',...
                    'is feasible using OPTIONS.TolCon of %1.6f.'],...
                    output.message, tolCon);
            end
            if do_display_P
                fprintf(1, output.message); end 
            
            go_on = false;
            return;
            
        end
                
        % If all variables are fixed, simply return
        if sum(lb(:)==ub(:)) == N0
            
            warning([mfilename ':bounds_overconstrain'],...
                    'Lower and upper bound are equal - nothing to do.');
            
            sol   = reshape(lb,size(x0));
            fval  = objFcn(sol);
            new_x = sol;
            
            exitflag = 2;
            if create_output
                output.iterations = 0;
                output.message = 'Lower and upper bound were set equal - nothing to do. ';
                incrementObjfuncCount();
            end
            
            do_display_P = do_display;
            do_display   = false;
            [output, exitflag] = finalize(sol, output, exitflag);
            
            if create_output && exitflag ~= -2
                output.message = sprintf(['%s\nFortunately, the solution ',...
                    'is feasible using OPTIONS.TolCon of %1.6f.'],...
                    output.message, tolCon);
            end            
            if do_display_P
                fprintf(1, output.message); end 
            
            go_on = false;
            return;
        end
             
    end % check_input
    
    
    % Counter management
    function incrementObjfuncCount(amt)
        if create_output  
            if nargin==0, amt = 1; end
            if ~have_nonlconFcn
                if ~isfield(output, 'funcCount')
                    output.funcCount = amt;
                else
                    output.funcCount = output.funcCount + amt;
                end
            else
                if ~isfield(output, 'ObjfuncCount')
                    output.ObjfuncCount = amt;
                else
                    output.ObjfuncCount = output.ObjfuncCount + amt;
                end
            end
        end
    end
    function incrementConfuncCount(amt)
        if create_output                  
            if nargin==0, amt = 1; end            
            if ~isfield(output, 'ConstrfuncCount')
                output.ConstrfuncCount = amt;                
            else
                output.ConstrfuncCount = output.ConstrfuncCount + amt;
            end
        end
    end
    
    
    % Evaluate objective function
    function varargout = objFcn(x)
        [varargout{1:nargout}] = feval(funfcn,...
                                       reshape(x,size(new_x)),...
                                       varargin{:});        
    end
    
    
    % Evaluate non-linear constraint function
    function [c,ceq, grad_c,grad_ceq] = conFcn(x)
        
        c        = [];
        ceq      = [];
        grad_c   = [];   
        grad_ceq = [];
        
        x = reshape(x,size(new_x));
        
        if have_nonlconFcn
            if nonlconFcn_in_objFcn
                if grad_nonlcon_from_nonlconFcn
                    if grad_obj_from_objFcn
                        [~,~, c, ceq, grad_c, grad_ceq] = objFcn(x);
                    else
                        [~, c, ceq, grad_c, grad_ceq] = objFcn(x);
                    end
                else
                    if grad_obj_from_objFcn
                        [~,~, c, ceq] = objFcn(x);
                    else
                        [~, c, ceq] = objFcn(x);
                    end
                end
            else
                if grad_nonlcon_from_nonlconFcn
                    [c, ceq, grad_c, grad_ceq] = feval(nonlcon, x, varargin{:});        
                else
                    [c, ceq] = feval(nonlcon, x, varargin{:});        
                end
            end
        end        
    end
    
    
    % Create transformed variable X to conform to upper and lower bounds
    function Z = X(x)
        
        % Initialize 
        if (size(x,2) == 1)
            Z = Nzero;    rep = 1;
        else
            Z = Np1zero;  rep = None;
        end  
        
        % First insert fixed values...
        y = x0one(:, rep);
        y( fix_var,:) = lb(fix_var,rep);
        y(~fix_var,:) = x;
        x = y;  
        
        % ...and transform.
        Z(lb_only, :) = lb(lb_only, rep) + x(lb_only, :).^2;
        Z(ub_only, :) = ub(ub_only, rep) - x(ub_only, :).^2;
        Z(fix_var, :) = lb(fix_var, rep);
        Z(unconst, :) = x(unconst, :);
        Z(lb_ub, :)   = lb(lb_ub, rep) + (ub(lb_ub, rep)-lb(lb_ub, rep)) .* ...
            (sin(x(lb_ub, :)) + 1)/2;        
    end % X
    
    
    % Derivatives of transformed X 
    function grad_Z = gradX(x, grad_x)
                
        % ...and compute gradient        
        grad_Z                       = grad_x(~fix_var);
        grad_Z(lb_only(~fix_var), :) = +2*grad_x(lb_only(~fix_var), :).*x(lb_only(~fix_var), :);
        grad_Z(ub_only(~fix_var), :) = -2*grad_x(ub_only(~fix_var), :).*x(ub_only(~fix_var), :);
        grad_Z(unconst(~fix_var), :) = grad_x(unconst(~fix_var), :);
        grad_Z(lb_ub(~fix_var), :)   = grad_x(lb_ub(~fix_var),:).*(ub(lb_ub(~fix_var),:)-lb(lb_ub(~fix_var),:)).*cos(x(lb_ub(~fix_var),:))/2;
        
    end % grad_Z
    
    
    % Create penalized function. Penalize with linear penalty function if 
    % violation is severe, otherwise, use exponential penalty. If the
    % 'strict' option has been set, check the constraints, and return INF
    % if any of them are violated.
    function [P_fval, grad_val] = funfcnP(x)
% FIXME: FMINSEARCH() or FMINLBFGS() see this as ONE function
% evaluation. However, multiple evaluations of both objective and nonlinear
% constraint functions may take place
        
        % Initialize function value
        if (size(x,2) == 1)
            P_fval = 0;
        else
            P_fval = None.'-1;
        end
                        
        % Initialize x_new array
        x_new = new_x;
        
        % Initialize gradient when needed
        if grad_obj_from_objFcn
            grad_val = zeros(size(x)); end
        
        % Evaluate every column in x
        for ii = 1:size(x,2)
            
            % Reshape x, so it has the same size as the given x0
            x_new(:) = x(:,ii);
            
            % Initialize
                           obj_gradient = 0;
            c      = [];   ceq          = []; 
            grad_c = [];   grad_ceq     = [];
            
            % Evaluate the objective function, taking care that also 
            % a gradient and/or contstraint function may be supplied
            if grad_obj_from_objFcn                
                if ~nonlconFcn_in_objFcn
                    [obj_fval, obj_gradient] = objFcn(x_new);
                    
                else
                    if grad_nonlcon_from_nonlconFcn
                          arg_out = cell(1, nonlconFcn_in_objFcn+1);
                    else, arg_out = cell(1, nonlconFcn_in_objFcn+3);
                    end
                    
                    [arg_out{:}] = objFcn(x_new);                  
                    obj_fval     = arg_out{1};
                    obj_gradient = arg_out{2};
                    c            = arg_out{nonlconFcn_in_objFcn+0};
                    ceq          = arg_out{nonlconFcn_in_objFcn+1};
                    
                    if grad_nonlcon_from_nonlconFcn
                        grad_c   = arg_out{nonlconFcn_in_objFcn+2};
                        grad_ceq = arg_out{nonlconFcn_in_objFcn+3};
                    else
                        grad_c   = ''; % use strings to distinguish them later on
                        grad_ceq = '';
                    end
                    
                end
                
                objFcn_fevals = 1;
                
            else
% FIXME: grad_nonlcon_from_nonlconFcn?
                if ~nonlconFcn_in_objFcn
                    obj_fval = objFcn(x_new);                  
                                        
                else                    
                    arg_out = cell(1, nonlconFcn_in_objFcn+1);
                    [arg_out{:}] = objFcn(x_new);                  
                    obj_fval     = arg_out{1};
                    c            = arg_out{nonlconFcn_in_objFcn+0};
                    ceq          = arg_out{nonlconFcn_in_objFcn+1};
                    grad_c       = ''; 
                    grad_ceq     = ''; % use strings to distinguish them later on
                    
                end
                
                objFcn_fevals = 1;
                
                if need_grad
                    [objFevals, obj_gradient] = computeJacobian(funfcn, x_new, obj_fval);
                    objFcn_fevals = objFcn_fevals + objFevals;
                end
                
            end  
            
            % Keep track of function evaluations
            incrementObjfuncCount(objFcn_fevals);
                                    
            % Make global copy
            UPfval = obj_fval; 
            
            % Initially, we are optimistic
            linear_eq_Penalty   = 0;    linear_ineq_Penalty_grad = 0;        
            linear_ineq_Penalty = 0;    linear_eq_Penalty_grad   = 0;
            nonlin_eq_Penalty   = 0;    nonlin_eq_Penalty_grad   = 0;            
            nonlin_ineq_Penalty = 0;    nonlin_ineq_Penalty_grad = 0;  
            
            % Penalize the linear equality constraint violation 
            % required: Aeq*x = beq   
            if have_lineqconFcn
                
                lin_eq = Aeq*x_new - beq;
                sumlin_eq = sum(abs(lin_eq(:)));
                
% FIXME: column sum is correct, but does not take into accuont non-violated
% constraints. We'll have to re-compute it 
                
                if strcmpi(strictness, 'All') && any(abs(lin_eq) > 0)
                    P_fval = inf; grad_val = inf; return; end
% FIXME: this really only works with fminsearch; fminlbfgs does not know
% how to handle this...
                
                % compute penalties
                linear_eq_Penalty = Penalize(sumlin_eq);
                
                % Also compute derivatives 
                % (NOTE: since the sum of the ABSOLUTE values is used
                % here, the signs are important!)
                if grad_obj_from_objFcn && linear_eq_Penalty ~= 0
                    linear_eq_Penalty_grad = ...
                        Penalize_grad(sign(lin_eq).*sumAeqT, sumlin_eq); 
                end
            end
                                                        
            % Penalize the linear inequality constraint violation 
            % required: A*x <= b
            if have_linineqconFcn
                
                lin_ineq = A*x_new - b;                     
                lin_ineq(lin_ineq <= 0) = 0;
                
% FIXME: column sum is correct, but does not take into accuont non-violated
% constraints. We'll have to re-compute it 
                
                sumlin_ineq      = sum(lin_ineq(:)); 
                sumlin_ineq_grad = sumAT;
                
                if strcmpi(strictness, 'All') && any(lin_ineq > 0)
                    P_fval = inf; grad_val = inf; return; end
% FIXME: this really only works with fminsearch; fminlbfgs does not know
% how to handle this...
                
                % Compute penalties
                linear_ineq_Penalty = Penalize(sumlin_ineq);
                
                % Also compute derivatives
                if grad_obj_from_objFcn && linear_ineq_Penalty ~= 0
                    linear_ineq_Penalty_grad = Penalize_grad(sumlin_ineq_grad, sumlin_ineq); end
                
            end
            
            % Penalize the non-linear constraint violations
            % required: ceq = 0 and c <= 0
            if have_nonlconFcn 
                
                [c, ceq, grad_c, grad_ceq] = conFcn(x_new);
                
                % Central-difference derivatives are computed later;
                % the strictness setting might make computing it here
                % unneccecary
                
                % Initialize as characters, to distinguish them later on;
                % derivatives may be empty, inf, or NaN as returned from [nonlcon]
                if isempty(grad_c)  , grad_c   = ''; end
                if isempty(grad_ceq), grad_ceq = ''; end
                
                %{
                if ~nonlconFcn_in_objFcn
                    
                    % Initialize as characters, to distinguish them later on;
                    % derivatives may be empty, inf, or NaN as returned from [nonlcon]
                    grad_c   = '';
                    grad_ceq = '';
                    
                    % Gradients are given explicitly by [nonlcon]
                    if grad_nonlcon_from_nonlconFcn
                        [c, ceq, grad_c, grad_ceq] = conFcn(x_new);
                        % The gradients are not given by [nonlcon]; they have to
                        % be computed by central differences
                    else
                        [c, ceq] = conFcn(x_new);
                        % Central-difference derivatives are computed later;
                        % the strictness setting might make computing it here
                        % unneccecary
                    end
                    
                    % Keep track of number of evaluations made
                    incrementConfuncCount();
                    
                else
                    % TODO
                end
                %}
                incrementConfuncCount();                
                
            end
            
            % Force grad_c] and [grad_ceq] to be of proper size
            if ~isempty(grad_c)
                grad_c = reshape(grad_c(:), numel(c), numel(x_new)); end
            if ~isempty(grad_ceq)
                grad_ceq = reshape(grad_ceq(:), numel(ceq), numel(x_new)); end

            % Process non-linear inequality constraints
            if ~isempty(c)                
                c = c(:); 
                
                % check for strictness setting
                if any(strcmpi(strictness, {'Bounds' 'All'})) &&...
                        any(c > ~superstrict*tolCon)
                    P_fval = inf; grad_val = inf; return
%FIXME: this only makes sense for FMINSEARCH
                end  
                
                % sum the violated constraints
                violated_c = c > tolCon;
                sumc = sum(c(violated_c));
                
                % compute penalty
                nonlin_ineq_Penalty = Penalize(sumc);
            end

            % Process non-linear equality constraints
            if ~isempty(ceq)
                % Use the absolute values, but save the signs for the
                % derivatives
                signceq = repmat(sign(ceq), 1,numel(x_new)); 
                ceq = abs(ceq(:));
                
                % Check for strictness setting
                if (strcmpi(strictness, 'Bounds') || ...
                        strcmpi(strictness, 'All')) &&...
                        any(ceq >= ~superstrict*tolCon)
                    P_fval = inf; grad_val = inf; return
%FIXME: this only makes sense for FMINSEARCH
                end 
                
                % Sum the violated constraints
                violated_ceq = (ceq >= tolCon); 
                sumceq = sum(ceq(violated_ceq));
                
                % Compute penalty 
                nonlin_eq_Penalty = Penalize(sumceq);
            end

            % Compute derivatives with central-differences of non-linear constraints
            if grad_obj_from_objFcn && ischar(grad_c) && ischar(grad_ceq)                
                [conFcn_fevals, grad_c, grad_ceq] = computeJacobian(nonlcon, x_new, c, ceq);
                incrementConfuncCount(conFcn_fevals);
            end
                        
            % Add derivatives of non-linear equality constraint function
            if grad_obj_from_objFcn && ~isempty(c)

                % First, remove those that satisfy the constraints
                grad_c = grad_c(violated_c, :);
                % Compute derivatives of penalty functions
                if ~isempty(grad_c)
                    nonlin_ineq_Penalty_grad = ...
                        Penalize_grad(sum(grad_c,1), sumc);
                end
            end

            % Add derivatives of non-linear equality constraint function
            if grad_obj_from_objFcn && ~isempty(ceq)

                % First, remove those that satisfy the constraints
                grad_ceq = grad_ceq(violated_ceq, :);
                % Compute derivatives of penalty functions
                % (NOTE: since the sum of the ABSOLUTE values is used
                % here, the signs are important!)
                if ~isempty(grad_ceq)
                    nonlin_eq_Penalty_grad = ...
                        Penalize_grad(sum(signceq.*grad_ceq,1), sumceq); 
                end
            end                
            
            % Return penalized function value
            P_fval(ii) = obj_fval + linear_eq_Penalty + linear_ineq_Penalty + ...
                nonlin_eq_Penalty + nonlin_ineq_Penalty; %#ok MLINT is wrong here...
            
             % Return penalized derivatives of constraints             
            grad_val(:, ii) = obj_gradient(:) + linear_eq_Penalty_grad(:) + ...
                linear_ineq_Penalty_grad(:) + nonlin_eq_Penalty_grad(:) + nonlin_ineq_Penalty_grad(:);
            
        end
        
        % Compute deserved penalties
        function fP = Penalize(violation)
            
            if violation <= tolCon
                fP = 0; return; end
            
            % Scaling parameter
% FIXME: scaling does not apear to work very well with fminlbfgs
            if strcmpi(algorithm, 'fminsearch')
                scale = min(1e60, violation/tolCon);          
            else
                scale = 1; 
            end
            
            % Linear penalty to avoid overflow
            if scale*violation > 50
                fP = exp50*(1 + scale*violation) - 1;
                
             % Exponential penalty otherwise
            else
                fP = exp(scale*violation) - 1;
            end
            
        end % Penalize
        
        % Compute gradient of penalty function
        % (doubly-nested function)
        function grad_fP = Penalize_grad(dvdx, violation)
           
            if violation <= tolCon
                grad_fP = 0; return; end
            
            % Scaling parameter
% FIXME: scaling does not apear to work very well with fminlbfgs
            dsdx = 0;
            if strcmpi(algorithm, 'fminsearch')
                scale = min(1e60, violation/tolCon);
                if violation/tolCon < 1e60
                    dsdx = dvdx/tolCon; end
            else
                scale = 1;
            end
            
            % Derivative of linear penalty function
            if scale*violation > 50
                grad_fP = exp50*(dsdx*violation + scale*dvdx);
                
            % Derivative of exponential penalty function
            else
                grad_fP = (dsdx*violation + scale*dvdx)*exp(scale*violation);
            end
            
        end % Penalize_grad
        
    end % funfcnP
    
    
    % Define the transformed & penalized function    
    function varargout = funfcnT(x)
        
        % Compute transformed variable
        XT = X(x);
        
        % WITH gradient
        if grad_obj_from_objFcn            
            [varargout{1}, grad_val] = funfcnP(XT);
            % Transform gradient and output
            varargout{2} = gradX(x, grad_val);
            
        % WITHOUT gradient
        else
            varargout{1} = funfcnP(XT);
        end 
        
    end % funfcnT
    
    
    % Compute gradient/Jacobian with finite differences
    function [fevals, varargout] = computeJacobian(F, x, varargin)
                        
% FIXME: does this make sense?        
        perturb = min(max(diffMinChange, 1e-6), diffMaxChange); 
        fevals  = 0;
        narg    = numel(varargin);
        
        % initialize Jacobian
        J = cellfun(@(y)zeros(numel(y), numel(x)), ...
                    varargin, ...
                    'UniformOutput', false);
        
        % And compute it using selected method
        switch lower(finDiffType)
            
            case 'forward'
                dx_plus = cell(narg,1);
                for jj = 1:numel(x)
                    
                    x(jj) = x(jj) + perturb;
                    [dx_plus{1:narg}] = F(x);
                    x(jj) = x(jj) - perturb;
                                        
                    for kk = 1:narg
                        newGrad = (dx_plus{kk}-varargin{kk})/perturb;
                        if ~isempty(newGrad)
                            J{kk}(:,jj) = newGrad; end
                    end
                    
                    fevals = fevals + 1;
                    
                end
               
            case 'backward'                
                dx_minus = cell(narg,1);
                for jj = 1:numel(x)
                    
                    x(jj) = x(jj) - perturb;
                    [dx_minus{1:narg}] = F(x);
                    x(jj) = x(jj) + perturb;
                                        
                    for kk = 1:narg
                        newGrad = (varargin{kk}-dx_minus{kk})/perturb;
                        if ~isempty(newGrad)
                            J{kk}(:,jj) = newGrad; end
                    end
                    
                    fevals = fevals + 1;
                    
                end                
                
            case 'central'                
                dx_plus  = cell(narg,1);
                dx_minus = cell(narg,1);
                for jj = 1:numel(x)
                    
                    % Forward
                    x(jj) = x(jj) + perturb;
                    [dx_plus{1:narg}] = F(x);
                    
                    % Backward
                    x(jj) = x(jj) - 2*perturb;
                    [dx_minus{1:narg}] = F(x);
                    
                    % Reset x
                    x(jj) = x(jj) + perturb;
                    
                    % Insert new derivatives
                    for kk = 1:narg
                        newGrad = (dx_plus{kk}-dx_minus{kk})/2/perturb;
                        if ~isempty(newGrad)
                            J{kk}(:,jj) = newGrad; end
                    end
                    
                    fevals = fevals + 2;
                    
                end
                
            case 'adaptive'
                % TODO
                
        end
        
        % Assign all outputs
        varargout = J;
        
    end
   
    
    % Simple wrapper function for output and plot functions; 
    % these need to be evaluated with the UNtransformed variables
    function stop = OutputFcn_wrapper(x, optimvalues, state)
        % Transform x        
        x_new = new_x;  x_new(:) = X(x);
        % Unpenalized function value
        optimvalues.fval = UPfval;
        % Evaluate all output functions        
        stop = zeros(size(OutputFcn));
        for ii = 1:numel(OutputFcn)
            stop(ii) = feval(OutputFcn{ii}, x_new, optimvalues, state); end
        stop = any(stop);        
    end % OutputFcn_wrapper
    
    
    function stop = PlotFcn_wrapper(x, optimvalues, state)
        % Transform x        
        x_new = new_x;  x_new(:) = X(x);
        % Unpenalized function value
        optimvalues.fval = UPfval;
        % Evaluate all plot functions
        stop = zeros(size(PlotFcn));
        for ii = 1:numel(PlotFcn)
            stop(ii) = feval(PlotFcn{ii}, x_new, optimvalues, state); end
        stop = any(stop);
    end % PlotFcn_wrapper
    
    
    % Finalize the output
    function [output, exitflag] = finalize(x, output, exitflag)
        
        % reshape x so it has the same size as x0
        x_new = new_x; x_new(:) = x;
            
        % compute violations (needed in both display and output structure)
                  
        % initialiy we're optimistic
        is_violated   = false;
        max_violation = 0;
        
        % add proper [constrviolation] field
        if have_linineqconFcn
            Ax        = A*x_new;
            violated  = Ax >= b + tolCon;
            violation = Ax - b;
            violation(~violated) = 0; clear Ax
            output.constrviolation.lin_ineq{1} = violated;
            output.constrviolation.lin_ineq{2} = violation;
            is_violated   = is_violated || any(violated(:));
            max_violation = max(max_violation, max(violation));
            clear violation violated
        end
        if have_lineqconFcn
            Aeqx = Aeq*x_new;
            violated  = abs(Aeqx - beq) > tolCon;
            violation = Aeqx - beq;
            violation(~violated) = 0; clear Aeqx
            output.constrviolation.lin_eq{1} = violated;
            output.constrviolation.lin_eq{2} = violation;
            is_violated   = is_violated || any(abs(violated(:)));
            max_violation = max(max_violation, max(abs(violation)));
            clear violation violated
        end
        if have_nonlconFcn
            [c, ceq] = conFcn(x_new);
            incrementConfuncCount();                         
            if ~isempty(ceq)
                violated = abs(ceq) > tolCon;
                ceq(~violated) = 0;
                output.constrviolation.nonlin_eq{1} = violated;
                output.constrviolation.nonlin_eq{2} = ceq;
                is_violated = is_violated || any(violated(:));
                max_violation = max(max_violation, max(abs(ceq)));
                clear violation violated ceq
            end
            if ~isempty(c)
                violated = c > tolCon;
                c(~violated) = 0;
                output.constrviolation.nonlin_ineq{1} = violated;
                output.constrviolation.nonlin_ineq{2} = c;
                is_violated = is_violated || any(violated(:));
                max_violation = max(max_violation, max(c));
                clear violation violated c
            end
            clear c ceq
        end
                
        % Adjust output message
        if create_output && exitflag == -3 
            output.message = sprintf(...
                ' No finite function values encountered.\n');
        end
        if ~isfield(output, 'message'), output.message = ''; end % (safeguard)
        if is_violated
            exitflag = -2;
            message = sprintf(...
                [' Unfortunately, the solution is infeasible for the given value ',...
                'of OPTIONS.TolCon of %1.6e\n Maximum constraint violation: ',...
                '%1.6e'], tolCon, max_violation);
            clear max_violation
        else
            if exitflag >= 1, message = sprintf('\b\n and'); 
            else,             message = sprintf('\b\n but');
            end
            message = [message, sprintf([' all constraints are satisfied using ',...
                    'OPTIONS.TolCon of %1.6e.'], tolCon)];
        end
        
        % Display or update output structure
        if create_output            
            output.message = sprintf('%s\n%s\n', output.message, message); end
        if do_display && ~do_global_opt
            fprintf(1, '%s\n', message); end
        
        % Correct for output possibly wrongfully created above
        if ~create_output, output = []; end
                
    end % finalize  
    
    
    % Optimize global problem
    function [sol, fval, exitflag, output, grad, hess] = minimize_globally()
                
        % First perform error checks
        if isempty(ub) || isempty(lb) || any(isinf(lb)) || any(isinf(ub))
            error([mfilename ':lbub_undefined'],...
                  ['When optimizing globally ([x0] is empty), both [lb] and [ub] ',...
                  'must be non-empty and finite.'])
        end  
        
        grad = {};
        hess = {};
        have_deriv = strcmpi(algorithm, 'fminlbfgs');
        
        % Global minimum (for output function)
        glob_min = inf;
        
        % We can give the popsize in the options structure, 
        % or we use 25*(number of dimensions) individuals by default
        popsize = getoptimoptions('popsize', 25*numel(lb));
                   
        % Initialize population of random starting points
        population = repmat(lb,[1,1,popsize]) + ...
            rand(size(lb,1),size(lb,2), popsize).*repmat(ub-lb,[1,1,popsize]);
               
        % Get options, and reset maximum allowable function evaluations
        %maxiters   = getoptimoptions('MaxIter', 200*numel(lb));
        maxfuneval = getoptimoptions('MaxFunEvals', 1e4);
        MaxFunEval = floor( getoptimoptions('MaxFunEvals', 1e4) / popsize / 1.2);
        options    = setoptimoptions(options, 'MaxFunEvals', MaxFunEval);        
                
        % Create globalized wrapper for outputfunctions
        have_glob_OutputFcn = false;
        if ~isempty(options.OutputFcn)
            have_glob_OutputFcn = true;
            glob_OutputFcn = options.OutputFcn;
            options.OutputFcn = @glob_OutputFcn_wrapper;            
        end
        
        % First evaluate output function
        if have_glob_OutputFcn
            optimValues.iteration = 0;
            optimValues.x         = x0;
            optimValues.fval      = glob_min;
            optimValues.procedure = 'init';
            optimValues.funcCount = 0;
            if have_nonlconFcn % constrained problems
                optimValues.ConstrfuncCount = 1; end % one evaluation in check_input()
            glob_OutputFcn(x0, optimValues, 'init');
        end
        
        % Display header
        if do_display
            if do_extended_display
                fprintf(1, ['  Iter  evals    min f(x)    global min f(x)   max violation\n',....
                    '=====================================================================\n']);
            else
                fprintf(1, ['  Iter  evals    min f(x)    global min f(x)\n',....
                    '============================================\n']);
            end
        end
        
        % Slightly loosen options for global method, and
        % kill all display settings
        global_options = setoptimoptions(options, ...
                                         'TolX'   , 1e2 * options.TolX,...
                                         'TolFun' , 1e2 * options.TolFun,...
                                         'display', 'off');  
                                
        % Initialize loop
        best_fval = inf; iterations = 0; obj_evals = 0; new_x = population(:,:,1);
        sol = NaN(size(lb)); fval = inf; exitflag = []; output = struct; con_evals = 0;        

        % Loop through each individual, and use it as initial value
        for ii = 1:popsize
                      
            % Optimize current problem
            if have_deriv
                [sol_i, fval_i,...
                exitflag_i, output_i,...
                ~,~] = minimize(funfcn,...
                                population(:,:,ii), ...
                                A,b,...
                                Aeq,beq,...
                                lb,ub,...
                                nonlcon,...
                                global_options);                
            else
                [sol_i, fval_i,...
                exitflag_i, output_i] = minimize(funfcn,...
                                                  population(:,:,ii), ...
                                                  A,b,...
                                                  Aeq,beq,...
                                                  lb,ub,...
                                                  nonlcon,...
                                                  global_options);
            end
    
            % Add number of evaluations and iterations to total
            if ~have_nonlconFcn % unconstrained problems
                obj_evals = obj_evals + output_i.funcCount;
            else % constrained problems
                obj_evals = obj_evals + output_i.ObjfuncCount;
                con_evals = con_evals + output_i.ConstrfuncCount;
            end   
            iterations = iterations + output_i.iterations;
            
            % Keep track of the best solution found so far
            if fval_i < best_fval
                % output values                
                fval = fval_i;   exitflag = exitflag_i;
                sol  = sol_i;    output   = output_i;
                % and store the new best
                best_fval = fval_i;
            end
            
            % Reset output structure
            if create_output
                if ~have_nonlconFcn % unconstrained problems
                    output.funcCount = obj_evals;
                else % constrained problems
                    output.ObjfuncCount = obj_evals;
                    output.ConstrfuncCount = con_evals;
                end
                output.iterations = iterations;
            end
            
            % Display output so far
            if do_display
                % iter-detailed: include max. constraint violation
                if do_extended_display
                    % do a dummy finalization to get the maximum violation
                    output_j = finalize(sol, output, exitflag_i);
                    max_violation = 0;
                    if have_nonlconFcn
                        max_violation = max([max_violation
                            output_j.constrviolation.nonlin_ineq{2}
                            abs(output_j.constrviolation.nonlin_eq{2})]);
                    end
                    if have_lineqconFcn
                        max_violation = max([max_violation
                            abs(output_j.constrviolation.lin_eq{2})]);
                    end
                    if have_linineqconFcn
                        max_violation = max([max_violation
                            output_j.constrviolation.lin_ineq{2}]);
                    end
                    % print everything
                    fprintf(1, '%4.0d%8.0d%14.4e%15.4e%16.4e\n', ...
                        ii,obj_evals, fval_i,best_fval,max_violation);
                % iter: don't
                else
                    % just print everything
                    fprintf(1, '%4.0d%8.0d%14.4e%15.4e\n', ...
                        ii,obj_evals, fval_i,best_fval);
                end
            end

            % MaxIters & MaxEvals check. The output function may also have 
            % stopped the global optimization
            % TODO: iterations...Document this change
            if (exitflag_i == -1) || ...%(iterations >= maxiters) || ...
               (obj_evals >= maxfuneval)
                % finalize solution                
                [dummy, exitflag] = finalize(sol, output, exitflag);%#ok
                % and break (NOT return; otherwise the last evaluation of 
                % the outputfunction will be skipped)
                break;
            end
            
        end % for
   
        % final evaluate output function
        if have_glob_OutputFcn
            optimValues.iteration = iterations;            
            optimValues.procedure = 'optimization complete';  
            optimValues.fval      = fval;
            optimValues.x         = sol;
            optimValues.funcCount = obj_evals;
            if have_nonlconFcn % constrained problems
                optimValues.ConstrfuncCount = con_evals; end
            glob_OutputFcn(sol, optimValues, 'done');
        end
        
        % check for INF or NaN values. If there are any, finalize 
        % solution and return
        if ~isfinite(fval)    
            [output, exitflag] = finalize(sol, output, -3); 
            if nargout >= 5
                grad = NaN(numel(lb),1);
                hess = NaN(numel(lb));
            end
            return; 
        end
         
        % Reset max. number of function evaluations
        options.MaxFunEvals = maxfuneval - obj_evals;
        if have_nonlconFcn % correction for constrained problems
            options.MaxFunEvals = maxfuneval - obj_evals-con_evals; end
        
        % Make 100% sure the display if OFF
        options.Display = 'off';
        
        % Perform the final iteration on the best solution found
        % NOTE: minimize with the stricter options
        if have_deriv
            [sol, fval,...
            exitflag, output_i,...
            grad,hess] = minimize(funfcn, sol,...
                                  A,b,...
                                  Aeq,beq,...
                                  lb,ub,...
                                  nonlcon,...
                                  options);
        else
            [sol, fval,...
            exitflag, output_i] = minimize(funfcn, sol,...
                                            A,b,...
                                            Aeq,beq,...
                                            lb,ub,...
                                            nonlcon,...
                                            options);
        end
        % Adjust output        
        if create_output
            if ~have_nonlconFcn % unconstrained problems
                output.funcCount = output.funcCount + output_i.funcCount;
            else % constrained problems
                output.ObjfuncCount    = output.ObjfuncCount + output_i.ObjfuncCount;
                output.ConstrfuncCount = output.ConstrfuncCount + output_i.ConstrfuncCount;
            end
            output.iterations = output.iterations + output_i.iterations;
        end
        
        % Get the final display right
        if do_display
            fprintf(1, output_i.message); end

        % Create temporary message to get the display right
        output.message = output_i.message;
        
        % Globalized wrapper for output functions
        function stop = glob_OutputFcn_wrapper(x, optimvalues, state) 
            % only evaluate if the current function value is better than
            % the best thus far found. Also evaluate on on first and last call
            stop = false;
            if (optimvalues.fval <= glob_min) &&...
                    ~any(strcmpi(state, {'done'; 'init'}))
                glob_min = optimvalues.fval;
                stop = glob_OutputFcn(x, optimvalues, state);            
            end
        end
        
    end % minimize_globally
    
    
    % Safe getter for (customized) options structure.
    % NOTE: this is a necessary workaround, because optimget() does not
    % handle non-standard fields.
    function value = getoptimoptions(parameter, defaultValue)
        
        if ~ischar(parameter)
            error('getoptimoptions:invalid_parameter',...
                'Expected parameter of type ''char'', got ''%s''.', class(parameter));
        end
        
        value = [];
        if isfield(options, parameter) && ~isempty(options.(parameter))
            value = options.(parameter);
        elseif nargin == 2
            value = defaultValue;
        end        
    end
    
    
end % function
