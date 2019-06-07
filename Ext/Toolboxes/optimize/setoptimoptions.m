% SETOPTIMOPTIONS                   Optimization options for minimize
%
% Usage:
%
%     options = setoptimoptions('param1',value1, 'param2',value2, ...)
%
% Options used by minimize() are:
%
%   Same as in <a href="matlab:doc('optimset')">optimset</a>:
%   - TolFun, TolX, OutputFcn, PlotFcn, MaxIter, MaxFunEvals, Display.
%
%   Specific to <a href="matlab:doc('minimize')">minimize</a>:
%       - TolCon
%           Tolerance used on any constraint. This means that minimize()
%           considers the constraints are only violated when they exceed
%           this amount of violation; otherwise, the constraints are
%           considered met. Defaults to 1e-8.
%
%       - GradObj
%           Specifies whether the objcetive function returns gradient
%           information as its second output argument. Valid values are
%           'on' and 'off' (the default). In case this option is set to
%           'off', gradient information is computed via finite
%           differences.
%
%       - GradConstr
%           Specifies whether the non-linear constraint function returns
%           Jacobian information as its third and fourth output arguments.
%           Valid values are 'on' and 'off' (the default). In case this
%           option is 'off', Jacobian information is computed via finite
%           differences.
%
%       - FinDiffType
%           Type of finite differences to use. Valid values are 'forward'
%           (the default), 'backward', and 'central'. Central differences
%           provide the best accuracy, but require twice as many function
%           evaluations.
%
%       - DiffMaxChange
%           Maximum change in the objective/constraint function to use when
%           computing gradient/Jacobian information with finite
%           differences. Defaults to 1e-1.
%
%       - DiffMinChange
%           Minimum change in the objective/constraint function to use when
%           computing gradient/Jacobian information with finite
%           differences. Defaults to 1e-8.
%
%       - AlwaysHonorConstraints       {'none'} 'bounds' 'all'
%           By default, minimize() will assume the objective (and
%           constraint) function(s) can be evaluated at ANY point in
%           RN-space; the initial estimate does not have to lie in the
%           feasible region, and intermediate solutions are also allowed
%           to step outside this area. this is equal to setting this option
%           to 'none'.
%
%           If the objective function cannot be evaluated outside the
%           feasible region, set this argument to 'bounds' (bound
%           constraints will never be broken) or 'all' (also the linear
%           constraints will never be broken). Note that the non-linear
%           constraints will remain satisfied within options.TolCon.
%
%           When using 'Bounds' or 'All', the initial estimate [x0]
%           MUST be feasible. If it is not feasible, an error is produced
%           before the objective function is ever evaluated.
%
%       - Algorithm
%           By default, this is set to MATLAB's own derivative-free
%           Nelder-Mead algorithm, implemented in <a href="matlab:doc('fminsearch')">fminsearch</a>.
%           minimize() supported another algotithm, <a
%           href="matlab:doc('fminlbfgs')">fminlbfgs</a>, a limited-memory,
%           Broyden/Fletcher/Goldfarb/Shanno optimizer. Use this algorithm
%           when your objective function has many free variables, e.g., [x]
%           is large.
%
%       - popsize
%           Used by the global optimization routine. This is the number of
%           randomized initial values to use, and thus the number of times
%           to repeat the call to minimize(). Defaults to 20× the number of
%           elements in [x0].
%
%   Specific to <a href="matlab:doc('fminlbfgs')">fminlbfgs</a>:
%       - GoalsExactAchieve
%           If set to 0, a line search method is used which uses a few
%           function calls to do a good line search. When set to 1 a normal
%           line search method with Wolfe conditions is used (default).
%
%       - HessUpdate
%           If set to 'bfgs', Broydenï¿½Fletcherï¿½Goldfarbï¿½Shanno
%			optimization is used (default), when the number of unknowns is
%			larger then 3000 the function will switch to Limited memory BFGS,
%			or if you set it to 'lbfgs'. When set to 'steepdesc', steepest
%			decent optimization is used.
%
%       - StoreN
%           Number of itterations used to approximate the Hessian,
%			in L-BFGS, 20 is default. A lower value may work better with
%			non smooth functions, because than the Hessian is only valid for
%			a specific position. A higher value is recommend with quadratic
%           equations.
%
%       - rho
%           Wolfe condition on gradient (c1 on wikipedia), default 0.01.
%
%       - sigma
%           Wolfe condition on gradient (c2 on wikipedia), default 0.9.
%
%       - tau1
%           Bracket expansion if stepsize becomes larger, default 3.
%
%       - tau2
%           Left bracket reduction used in section phase, default 0.1.
%
%       - tau3
%           Right bracket reduction used in section phase, default 0.5.
%
% See also optimset, optimget.


% Please report bugs and inquiries to:
%
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com
% Licence    : 3-clause BSD (See License.txt)


% Changelog
%{
2018/June
- CHANGED: removed affilication info & updated license

2014/July/07 (Rody Oldenhuis)
- FIXED: loop range for toolbox options was set to the wrong options list.
         Problem reported by pag (https://github.com/pag). Thanks!

2014/March/13 (Rody Oldenhuis)
- Finally did the docs!

2014/February/19 (Rody Oldenhuis)
- initial version
%}

% TODO
%{
check in what MATLAB version inputParser() was introduced
%}
function options = setoptimoptions(varargin)

    % If you find this work useful, please consider a donation:
    % https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=6G3S5UYM7HJ3N

    %% Initialize

    % NOTE: a lot of this can be done with inputParser(). It has been
    % implemented like this to support older versions of MATLAB

    % Check if we have the "advanced" optimset/optimget
    persistent haveOptimToolbox
    if isempty(haveOptimToolbox)
        haveOptimToolbox = ~isempty(ver('optim')); end

    oldOptions = [];
    if nargin >= 1 && isstruct(varargin{1})
        oldOptions = varargin{1};
        parameters = varargin(2:2:end);
        values     = varargin(3:2:end);

    else
        % Basic check on inputs
        if mod(nargin,2)~=0
            error('setoptimoptions:pvpairs_expected',...
                'setoptimoptions expects parameter-value pairs.');
        end

        parameters = varargin(1:2:end);
        values     = varargin(2:2:end);
    end

    %% Delegate part of the work to optimset()

    % TODO: GradConstr means something else in fminlbfgs() than in minimize()

    % Custom options
    customOptions = {
        %{
        minimize() native
        %}
        'AlwaysHonorConstraints'        , 'none'    % "strictness"
        'Algorithm'                     , 'fminsearch'
        'ConstraintsInObjectiveFunction', false
        'FinDiffType'                   , 'forward' % forward, backward, central, adaptive,
        'popsize'                       , []
        %{
        Specific to fminlbfgs()
        %}
        'GoalsExactAchieve'             , 1         % Normal line search with Wolfe conditions
        'HessUpdate'                    , 'bfgs'    % bfgs, lbfgs, steepdesc
        'StoreN'                        , 20        % Number of itterations used to approximate the Hessian
        'rho'                           , 0.01      % Wolfe condition on gradient (c1 on wikipedia)
        'sigma'                         , 0.9       % Wolfe condition on gradient (c2 on wikipedia)
        'tau1'                          , 3         % Bracket expansion if stepsize becomes larger
        'tau2'                          , 0.1       % Left bracket reduction used in section phase
        'tau3'                          , 0.5       % Right bracket reduction used in section phase
    };


    customParameter = false(size(parameters));
    for ii = 1:size(customOptions,1)
        customParameter = customParameter | strcmpi(parameters, customOptions{ii,1}); end

    % Options only in the optimization toolbox version of optimset
    if ~haveOptimToolbox
        toolboxOptions = {
            'TolCon'         , 1e-8
            'GradObj'        , []
            'GradConstr'     , []
            'DerivativeCheck', 'off'
            'FinDiffType'    , 'forward'
            'DiffMaxChange'  , 1e-1
            'DiffMinChange'  , 1e-8
            };

        for ii = 1:size(toolboxOptions,1)
            customParameter = customParameter | strcmpi(parameters, toolboxOptions{ii,1}); end
    end

    % Delegate all non-custom parameters to optimset
    delegate = [parameters(~customParameter); values(~customParameter)];
    options  = optimset(delegate{:});

    % Initialize all custom parameters to their default values
    for ii = 1:size(customOptions,1)
        options.(customOptions{ii,1}) = customOptions{ii,2}; end
    if ~haveOptimToolbox
        for ii = 1:size(toolboxOptions,1)
            options.(toolboxOptions{ii,1}) = toolboxOptions{ii,2}; end
    end

    % Merge any old options with new ones
    if ~isempty(oldOptions)

        fOld = fieldnames(oldOptions);

        % Custom options
        for ii = 1:numel(customOptions(:,1))
            if ~any(strcmpi(parameters, customOptions{ii,1})) && isfield(oldOptions, customOptions{ii,1}) && ~isempty(oldOptions.(customOptions{ii,1}))
                options.(customOptions{ii,1}) = oldOptions.(customOptions{ii,1});
                fOld(strcmpi(fOld,customOptions{ii,1})) = [];
            end
        end

        % Toolbox options
        if ~haveOptimToolbox
            for ii = 1:numel(toolboxOptions(:,1))
                if ~any(strcmpi(parameters, toolboxOptions{ii,1})) && isfield(oldOptions, toolboxOptions{ii,1}) && ~isempty(oldOptions.(toolboxOptions{ii,1}))
                    options.(toolboxOptions{ii,1}) = oldOptions.(toolboxOptions{ii,1});
                    fOld(strcmpi(fOld,toolboxOptions{ii,1})) = [];
                end
            end
        end

        % All other options
        if ~isempty(fOld)
            for ii = 1:numel(fOld)
                if ~any(strcmpi(parameters, fOld{ii})) && isfield(oldOptions, fOld{ii}) && ~isempty(oldOptions.(fOld{ii}))
                    options.(fOld{ii}) = oldOptions.(fOld{ii}); end
            end
        end

    end


    %% Parse the PV-pairs

    % The remainder
    parameters = parameters(customParameter);
    values     = values(customParameter);

    if ~isempty(parameters)
        for ii = 1:numel(parameters)

            parameter = parameters{ii};
            value     = values{ii};

            switch lower(parameter)


                %% Optimize native

                case 'algorithm'
                    if ~isValidString(value, {'fminsearch', 'fminlbfgs'}, parameter)
                        continue; end
                    options.Algorithm = value;

                case 'alwayshonorconstraints'
                    if ~isValidString(value, {'none', 'bounds', 'all'}, parameter)
                        continue; end
                    options.AlwaysHonorConstraints = value;

                case 'constraintsinobjectivefunction'
                    if ~isClass('logical', value)
                        continue; end
                    options.ConstraintsInObjectiveFunction = value;

                case 'findifftype'
                    if ~isValidString(value, {'forward', 'backward', 'central', 'adaptive'}, parameter)
                        continue; end
                    options.FinDiffType = value;

                case 'popsize'
                    if ~isClass('numeric', value)
                        continue; end
                    options.popsize = value;


                %% Fminlbfgs

                case 'goalsexactachieve'
                    if ~isClass('numeric', value)
                        continue; end
                    options.GoalsExactAchieve = isfinite(value(1)) && value(1)~=0;

                case 'hessupdate'
                    if ~isValidString(value, {'bfgs', 'lbfgs', 'steepdesc'}, parameter)
                        continue; end
                    options.HessUpdate = value;

                case 'storen'
                    if ~isClass('numeric', value)
                        continue; end
                    options.StoreN = abs(value(1));

                case 'rho'
                    if ~isClass('numeric', value)
                        continue; end
                    options.rho = value(1);

                case 'sigma'
                    if ~isClass('numeric', value)
                        continue; end
                    options.sigma = value(1);

                case 'tau1'
                    if ~isClass('numeric', value)
                        continue; end
                    options.tau1 = value(1);

                case 'tau2'
                    if ~isClass('numeric', value)
                        continue; end
                    options.tau2 = value(1);

                case 'tau3'
                    if ~isClass('numeric', value)
                        continue; end
                    options.tau3 = value(1);


                %% Optimization toolbox

                otherwise

                    if ~haveOptimToolbox
                        switch lower(parameter)
                            case 'tolcon'
                                if ~isClass('numeric', value)
                                    continue; end
                                options.TolCon = value(1);

                            case 'gradobj'
                                % TODO: support function handles?
                                % TODO: Support string functions?
                                if ~isValidString(value, {'on', 'off'}, parameter)
                                    continue; end
                                options.GradObj = value;

                            case 'gradconstr'
                                % TODO: support function handles?
                                % TODO: Support string functions?
                                if ~isValidString(value, {'on', 'off'}, parameter)
                                    continue; end
                                options.GradConstr = value;

                            case 'derivativecheck'
                                if ~isValidString(value, {'yes', 'no'}, parameter)
                                    continue; end
                                % TODO: valid strings
                                options.DerivativeCheck = value;

                            case 'findifftype'
                                if ~isValidString(value, {'central', 'forward', 'backward'}, parameter)
                                    continue; end
                                options.FinDiffType = value;

                            case 'diffmaxchange'
                                if ~isClass('numeric', value)
                                    continue; end
                                options.DiffMaxChange = value(1);

                            case 'diffminchange'
                                if ~isClass('numeric', value)
                                    continue; end
                                options.DiffMinChange = value(1);

                            otherwise
                                warning(...
                                    'setoptimoptions:unknown_option',...
                                    'Unknown option: ''%s''. Ignoring...', parameter);
                                continue
                        end
                    else
                        warning(...
                            'setoptimoptions:unknown_option',...
                            'Unknown option: ''%s''. Ignoring...', parameter);
                        continue
                    end

            end % switch
        end % for
    end % if

end % function


% Check if given string is one of the supported options. Generate
% appropriate warning messages if this is not the case.

% NOTE: this is virtually the same as validatestring(), but this version
% can also be used on older versions of MATLAB
function go_on = isValidString(received, validStrings, parameter)

    go_on = true;

    if ~isClass('char', received)
        go_on = false; return; end

    if ~any(strcmpi(received, validStrings))

        validStrings = [
            cellfun(@(x) ['''' x ''', '], validStrings(1:end-1), 'UniformOutput', false),...
            ['and ''' validStrings(end) '''.']
        ];
        warning(...
            'setoptimoptions:unsupported_value',...
            ['Unsupported option ''%s'' for parameter ''%s''.\n',...
            'Valid options are: ' validStrings{:} '\nUsing default...'], ...
            received, parameter);

        go_on = false;
    end

end


% Check if class of given value is correct. Generate appropriate warning
% message if this is not the case.

% NOTE: this is virtually the same as validateattributes(), but this
% version can also be used on older versions of MATLAB
function go_on = isClass(expected, received)

    go_on = true;

    % More than one type may be allowable:
    if ischar(expected)
        isCorrectClass = isa(received, expected);

    elseif iscell(expected)
        isCorrectClass = any(cellfun(@(x) isa(received, x), expected));

    else
        error(...
            'setoptimoptions:BUG',...
            'An internal error has occurred; please report this to the lead developer.');
    end

    % Generate warning if the type is not any of the supported ones
    if ~isempty(received) && ~isCorrectClass
        warning(...
            'setoptimoptions:invalid_value',...
            'Expected ''%s'', got ''%s''. Ignoring...', expected, class(received));
        go_on = false;
    end

end


% Check validity of file/function given as strnig
% NOTE: this functionality is deprecated, so a warning is ALWAYS issued.
function [go_on, value] = checkFunction(value)

    go_on = true;

    if ischar(value)

        warning(...
            'setoptimoptions:deprecated',...
            'Using string functions or file names of MATLAB functions is deprecated; please use function_handles.');

        if ~any(exist(value,'file')==[2 3 5 6])
            warning(...
                'setoptimoptions:invalid_script',...
                ['The file/function ''%s'' is invalid, or is not on the MATLAB search path.\n',...
                'Using default...'], value);
            go_on = false;
        end

        % Strip any extension, and convert to valid function handle
        value = str2func(regexprep(value, '\.\w*$', ''));
    end
end

