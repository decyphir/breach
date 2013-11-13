function InitSensi(Sys, P)
%INITSENSI Initializes Cvodes for sensitivity computation of system Sys
% w.r.t. varying parameters in the set P. Options can be changed in the
% field CVodesSensiOptions
% 

dims = P.dim;
N = Sys.DimX;
ix0 = dims(dims<=N); % Uncertains initial conditions

%NM: I comment the following line, uncomment if something bad happens
%ip = dims(dims>N); % Uncertains parameters

InitSystem(Sys);

Ns = numel(dims);

if isfield(Sys,'CVodesSensiOptions')
    FSAoptions = Sys.CVodesSensiOptions.FSAoptions;
    FSAoptions = CVodeSetFSAOptions(FSAoptions,'ParamList', dims,...
        'ParamScales',FSAoptions.ParamScales(dims));
    method = Sys.CVodesSensiOptions.method;
    
else
    AbsTol = Sys.CVodesOptions.AbsTol;
    if isscalar(AbsTol)
        pscales = AbsTol*1e6*ones(Ns,1);
    else
        if ~isempty(ix0)
            pscales(dims<=N) = AbsTol(dims<=N)*1e6;
        else
            pscales(dims>N) = 1;
        end
    end
    method = 'Simultaneous';
    FSAoptions = CVodeSetFSAOptions('SensErrControl', 'on',...
        'ParamField', 'p',...
        'ParamList', dims,...
        'ParamScales',pscales );
end

xS0 = zeros(N,Ns);

for ii=1:numel(ix0);
    xS0(dims(ii),ii) = 1;
end;


CVodeSensMalloc(Ns, method, xS0, FSAoptions);

end
