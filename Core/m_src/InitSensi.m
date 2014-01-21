function InitSensi(Sys, P)
%INITSENSI initializes Cvodes for sensitivity computation of system Sys
% w.r.t. varying parameters in the set P. Options can be changed in the
% field CVodesSensiOptions
% 
% Synopsis: InitSensi(Sys, P)
%

dims = sort(P.dim);
N = Sys.DimX;
ix0 = dims(dims<=N); % Uncertains initial conditions

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
        pscales = AbsTol*1e6*ones(1,Ns); %NM: guess that pscales is a row
    else
        pscales = ones(1,Ns);
        pscales(1:numel(ix0)) = AbsTol(ix0)*1e6;
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
