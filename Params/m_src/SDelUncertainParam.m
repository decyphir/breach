function P = SDelUncertainParam(P, ParamList, Exception, Verbose)
%SDELUNCERTAINPARAM removes uncertains parameters of ParamList in P
% excepted of those contained in Exception. If P will contain no uncertain
% parameters, nothing is done.
% 
% Synopsis: P = SDelUncertainParam(P, ParamList[, Exception[, Verbose]])
% 
% Inputs:
%   - P         : the parameter set to modify
%   - ParamList : indexes or names of parameters to remove from P (excepted
%                 if they are also in Exception). Indexes or names not
%                 valid are ignored.
%   - Exception : (Optional, default=[]) indexes or names of parameters not
%                 to remove. Indexes or names not valid are ignored.
%   - Verbose   : (Optional, default=1). If set to one, the function shows
%                 a warning if the uncertain parameters are not removed.
% 
% Output:
%   - P : the modified parameter set
% 
% Examples (Lorentz84):
%   CreateSystem
%   P = CreateParamSet(Sys,{'a','b'});
%   fprintf('%s\n',P.ParamList{P.dim}); % output: a and b
%   P = SDelUncertainParam(P,{'a','F','blah'});
%   fprintf('%s\n',P.ParamList{P.dim});  % output: only b
% 
%See also SAddUncertainParam CreateParamSet SetEpsi GetEpsi
%

% Check inputs

if(iscell(ParamList) || ischar(ParamList))
    ParamList = FindParam(P,ParamList);
end

if ~exist('Exception','var')
    Exception = [];
elseif(iscell(Exception) || ischar(Exception))
    Exception = FindParam(P,Exception);
end

if ~exist('Verbose','var')
    Verbose = 1;
end

% does not remove Exception

idx_remove = ~ismember(ParamList,Exception); % logical indexes of uncertain removable parameters
idx = ~ismember(P.dim,ParamList(idx_remove));
if any(idx)
    P.epsi = P.epsi(idx,:);
    P.dim = P.dim(idx);
elseif(Verbose>=1)
    warning('SDelUncertainParam:NoUncertainParam',...
            ['No uncertain parameters are remaining, the command has' ...
             ' been skipped\n']);
end

end
