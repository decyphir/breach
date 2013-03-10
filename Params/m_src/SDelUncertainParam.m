function P = SDelUncertainParam(P, ParamList, Exception, Verbose)
%SDELUNCERTAINPARAM removes uncertains parameters of ParamList in P
% excepted of those contained in Exception. If P will contain no uncertain
% parameters, nothing is done.
%
% Synopsis: P = SDelUncertainParam(P, ParamList [ , Exception, Verbose ])
%
% Input:
%   - P         : the parameter set to modify
%   - ParamList : indexes or names of parameters to remove from P (excepted
%                 if they are also in Exception). Indexes or names not
%                 valid are ignored.
%   - Exception : (Optional) Indexes or names of parameters not to remove.
%                 Indexes or names not valid are ignored. (default=[])
%   - Verbose   : If set to one, the function shows a warning if the
%                 uncertain parameters are not removed (default=1).
% Output:
%   - P : the modified parameter set
%
% Examples (Lorentz84):
%
%   CreateSystem;
%   P = CreateParamSet(Sys,{'a','b'});
%   fprintf('%s\n',P.ParamList{P.dim}); % output: a and b
%   P = SDelUncertainParam(P,{'a','F','blah'});
%   fprintf('%s\n',P.ParamList{P.dim});  % output: only b
%
%See also SAddUncertainParam, CreateParamSet
%

% Check inputs

if iscell(ParamList) || ischar(ParamList)
    ParamList = FindParam(P,ParamList);
end

if iscell(Exception) || ischar(Exception)
    Exception = FindParam(P,Exception);
end

if ~exist('Exception','var')
    Exception = [];
elseif iscell(Exception) || ischar(Exception)
    Exception = FindParam(P,Exception);
end

if ~exist('Verbose','var')
    Verbose = 1;
end

% does not remove Exception

idx_remove = setdiff(ParamList,Exception,'stable'); % indexes of uncertain parameters to remove
[ParamList, idx] = setdiff(P.dim,idx_remove,'stable'); % indexes of uncertain parameters to keep
if ~isempty(ParamList)
    P.epsi = P.epsi(idx,:);
    P.dim = ParamList;
elseif(Verbose>=1)
    warning('SDelUncertainParam:NoUncertainParam',...
            ['No uncertain parameters are remaining, the command has' ...
             ' been skipped\n']);
end

% for ii=ParamList
%     if isempty(find(Exception==ii,1))
%         ParamList = find(P.dim~=ii);
%         P.dim = P.dim(ParamList);
%         P.epsi = P.epsi(ParamList,:);
%     end
% end

end
