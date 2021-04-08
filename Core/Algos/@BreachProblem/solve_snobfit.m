function res = solve_snobfit(this, varargin)
% SNOBFIT
% Install latest SNOBFIT, but MINQ for Matlab 5
% The following is mostly copied from
% snobfit/snobtest.m (some things changed)
file = 'snobfit_data';
fcn = 'snobfit_wrapper';
fac = 0;        % factor for multiplicative perturbation of the data
ncall = this.max_obj_eval;   % limit on the number of function calls

% Extract lower and upper bounds from this.domains
allDomains = {this.domains.domain};
lowerBounds = cellfun(@(x)x(1), allDomains);
upperBounds = cellfun(@(x)x(2), allDomains);
u = lowerBounds';
v = upperBounds';

startSample = this.solver_options.start_sample;
startFunctionValues = this.solver_options.start_function_values;
fglob = -0.01;
n = length(u);  % dimension of the problem
% the following are meaningful default values
npoint = 1;   % number of random start points to be generated
nreq = n+6;     % no. of points to be generated in each call to SNOBFIT
if isempty(startSample)
    % No startSample given
    startSample = testronGetNewSample([this.lb this.ub]);
else
    % Start sample is already given in the variable
    % startSample as an input to this function
end
x = startSample';
% starting points in [u,v]
dx = (v-u)'*1.e-5; % resolution vector
p = 0.5;        % probability of generating a point of class 4
prt = 0;        % print level
% prt = 0 prints ncall, xbest and fbest if xbest has
%         changed
% prt = 1 in addition prints the points suggested by
%         SNOBFIT, their model function values and
%         classes after each call to SNOBFIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of data to be adapted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(startFunctionValues)
    % No start function values provided
    for j=1:npoint
        functionEvalPlusNoise = feval(fcn,x(j,:),this)+fac*randn;
        if numel(functionEvalPlusNoise) > 1
            % We are evaluating several specs at once
            % Take the minimum of the robustness values
            % (i.e. use max semantics for the conjunction
            % of specs)
            functionEvalPlusNoise = min(functionEvalPlusNoise);
        end
        f(j,:) = [functionEvalPlusNoise max(sqrt(eps),3*fac)];
        % computation of the function values (if necessary, with additive
        % noise)
    end
    ncall0 = npoint;   % function call counter
else
    % Start function values provided!
	
	% Replace Inf values in startFunctionValues by really high numbers
    % snobfit cannot handle Inf
    startFunctionValues(isinf(startFunctionValues)) = 10000*max(this.avgRobForNormalization);
	
    f = [startFunctionValues' repmat(max(sqrt(eps),3*fac), numel(startFunctionValues), 1)];
    disp(['Starting SNOBFIT with ' num2str(size(f, 1)) ' pre-calculated objective function values']);
    ncall0 = 0;
    npoint = 0;
end
params = struct('bounds',{u,v},'nreq',nreq,'p',p); % input structure
% repeated calls to Snobfit
bestForPrint = Inf;
printFlag = 1;

while ncall0 < ncall % repeat till ncall function values are reached
    % (if the stopping criterion is not fulfilled first)
    if ncall0 == npoint  % initial call
        disp(['Initial robustness value: ' num2str(f(1))]);
        
        [request,xbest,fbest] = snobfit(file,x,f,params,dx);
        %ncall0,xbest,fbest;
    else                 % continuation call
        [request,xbest,fbest] = snobfit(file,x,f,params);
    end
    if prt>0, request, end
    clear x
    clear f
    
    if isfield(this.BrSys.Sys, 'use_parallel') && this.BrSys.Sys.use_parallel
        % Parallel computation
        this.use_parallel = 1;
        x = request(:, 1:n);
        fcn_handle = @snobfit_wrapper;
        % We utilize the parallel computation support of
        % BreachProblem.objective_wrapper
        allFunctionValues = snobfit_wrapper(x, this);
        for j = 1:size(allFunctionValues, 2)
            functionEvalPlusNoise = allFunctionValues(:, j) + fac*randn;
            if numel(functionEvalPlusNoise) > 1
                % We are evaluating several specs at once
                % Take the minimum of the robustness values
                % (i.e. use max semantics for the conjunction
                % of specs)
                functionEvalPlusNoise = min(functionEvalPlusNoise);
            end
            
            if isinf(functionEvalPlusNoise)
                % Infinite function value - snobfit can't handle
                % Replace infinite value by 10000 times the meximum
                % average robustness from normalization
                functionEvalPlusNoise = ...
                    10000*max(this.avgRobForNormalization);
            end
            
            f(j,:) = [functionEvalPlusNoise max(sqrt(eps),3*fac)];
            
            totalCounter = ncall0 + j - 1;
            if f(j,1) < bestForPrint
                bestForPrint = f(j,1);
                if ncall0 > npoint
                    disp([num2str(totalCounter) ': NEW BEST: ' num2str(bestForPrint)]);
                end
                if bestForPrint < 0
                    disp(['FALSIFIED at sample ' num2str(totalCounter) '!']);
                    printFlag = 0;
                end
            elseif mod(totalCounter, this.freq_update)==0 && printFlag && ~this.stopping
                fprintf([num2str(totalCounter) ': Rob: ' num2str(f(j,1)) '\t\tBEST:' num2str(bestForPrint) '\n']);
            end
        end
%         for j = 1:size(request, 1)
%             x(j,:) = parRequest(j, :);
%             parResults(j) = parfeval(fcn_handle, 1, x(j, :), this);
%         end
%         
%         for j2 = 1:size(request, 1)
%             [~, functionEvalPlusNoise] = fetchNext(parResults);
%             functionEvalPlusNoise = functionEvalPlusNoise + fac*randn;
%             
%             % Update obj_log and nb_obj_eval
%             this.obj_log = [this.obj_log functionEvalPlusNoise];
%             this.nb_obj_eval = this.nb_obj_eval + 1;
%             
%             if numel(functionEvalPlusNoise) > 1
%                 % We are evaluating several specs at once
%                 % Take the minimum of the robustness values
%                 % (i.e. use max semantics for the conjunction
%                 % of specs)
%                 functionEvalPlusNoise = min(functionEvalPlusNoise);
%             end
%             f(j2,:) = [functionEvalPlusNoise max(sqrt(eps),3*fac)];
%         end

    else
        % Serial computation
        for j=1:size(request,1)
            x(j,:) = request(j,1:n);
            functionEvalPlusNoise = feval(fcn,x(j,:),this)+fac*randn;
            if numel(functionEvalPlusNoise) > 1
                % We are evaluating several specs at once
                % Take the minimum of the robustness values
                % (i.e. use max semantics for the conjunction
                % of specs)
                functionEvalPlusNoise = min(functionEvalPlusNoise);
            end
            
            if isinf(functionEvalPlusNoise)
                % Infinite function value - snobfit can't handle
                % Replace infinite value by 10000 times the meximum
                % average robustness from normalization
                functionEvalPlusNoise = ...
                    10000*max(this.avgRobForNormalization);
            end
            
            f(j,:) = [functionEvalPlusNoise max(sqrt(eps),3*fac)];
            % computation of the (perturbed) function values at the suggested points
            
            % Display
            totalCounter = ncall0 + j - 1;
            if f(j,1) < bestForPrint
                bestForPrint = f(j,1);
                if ncall0 > npoint
                    disp([num2str(totalCounter) ': NEW BEST: ' num2str(bestForPrint)]);
                end
                if bestForPrint < 0
                    disp(['FALSIFIED at sample ' num2str(totalCounter) '!']);
                    printFlag = 0;
                end
            elseif mod(totalCounter, this.freq_update)==0 && printFlag && ~this.stopping
                fprintf([num2str(totalCounter) ': Rob: ' num2str(f(j,1)) '\t\tBEST:' num2str(bestForPrint) '\n']);
            end
        end
    end
    
    ncall0 = ncall0 + size(f,1); % update function call counter
    [fbestn,jbest] = min(f(:,1)); % best function value
    if fbestn < fbest
        
        fbest = fbestn;
        xbest = x(jbest,:);
        %ncall0,xbest,fbest % display current number of function values,
        % best point and function value if fbest has
        % changed
    end
    % check stopping criterion
    % if fglob == 0, stop if fbest < 1.e-5
    % otherwise, stop if (fbest-fglob)/abs(fglob) < 1.e-2
    if fbest < 0
        break
    end
    
    
end
%ncall0,xbest,fbest  % show number of function values, best point and
% function value
res = struct('bestRob',[],'bestSample',[],'nTests',[],'bestCost',[],'paramVal',[],'falsified',[],'time',[]);
res.bestSample = xbest;
res.bestRob = fbest;

end