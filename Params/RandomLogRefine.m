function [PRLog] = RandomLogRefine(P, N)
%RANDOMLOGREFINE Create a logarithmic random sampling of parameters. If P
%   contains many points, each of them is divided into N new points. All
%   ranges (ie : [value-epsi ; value+epsi] must be strictly positive. When
%   there is more than 16 ordre of magnitude (with Matlab R2012), it is
%   recommended to use CreateRandomLogParamSets or you risk to face error
%   due to null interval limit.
%
%   Syntax: PRLog = RandomLogSampling(P, N)
%
%   Inputs:
%
%    -  P         The parameter set to refine
%    -  N         Number of random generated points
%
%   Outputs:
%
%    -  PRLog    A random logarithmic sampling of N points
%
%	Example:
%
%   CreateSystem;
%   P = CreateParamSet(Sys);
%   PRLog = RandomLogRefine(P,10);
%
% SEE ALSO CREATERANDOMLOGPARAMSETS
%

% A AMELIORER :
%
% Utiliser la génération de nombre aléatoire type quasi monte-carlo au lieu
% du rand proposé par matlab.


nbParam = numel(P.dim);
params = cell(1,nbParam); %params est organisé dans l'ordre de P.dim
ranges = zeros(nbParam,2);
for i=1:nbParam
    params(i) = P.ParamList(P.dim(i));
end

% check for duplicates in params - probably not mandatory here
if(nbParam~=numel(unique(params)))
    %return ;
    error('RandomLogRefine:duplicateParam',['When refining, it must '...
        'not be more than one time a parameter']);
end

indicePts = 1; % l'indice de colonne auquel on ajoute les points et epsi calculé


PRLog = P;
PRLog.pts = zeros(size(P.pts,1),size(P.pts,2)*N);  %On découpe en N*le nombre de points dans P
PRLog.epsi = zeros(nbParam,size(P.pts,2)*N);

%we copy from P to PRLog all the parameter value which are not incertain (which don't change)
SD = setdiff(1:size(P.pts,1),P.dim); %index of not changing parameters
PRLog.pts(SD,:) = repmat(P.pts(SD,1),1,size(P.pts,2)*N);


for ii = 1:size(P.pts,2) % pour chacun des points de P, on découpe en N
    
    % on créé les ranges
    for i=1:nbParam
        ranges(i,1) = P.pts(P.dim(i),ii) - P.epsi(i,ii);
        ranges(i,2) = P.pts(P.dim(i),ii) + P.epsi(i,ii);
    end
    
    % we check if all range limits are strictly positive
    if(~isempty(find(sign(ranges)<=0,1)))
        %return;
        error('RandomLogSampling:rangeBound','Range limits must be strictly positive.');
    end
    
    ranges = log10(ranges);
    
    for i=1:nbParam
        %Pour chaque param, on défini une position aléatoire
        epsi = (ranges(i,2)-ranges(i,1))/2; %epsi de l'intervalle initial en log
        r = ranges(i,1) + 2*epsi*rand(1,N); %on tire au hasard sur l'échelle log
        epsi = epsi/(N^(1/nbParam)); %epsi des nouveaux intervalles en log
        
        % possibility 1 : define epsi, using the lower bound
        inf = r - epsi; %valeur minimal des nouveaux intervalles en log
        value = 10.^r; %valeur des nouveaux points en normal
        inf = 10.^inf; %valeur minimal des nouveaux intervalles en normal
        epsi = value - inf; %epsi des nouveaux intervalles en normal
        
        %possiblity 2 : define epsi using the higher bound
        %    sup = r + epsi; %valeur maximale des nouveaux intervalles en log
        %    value = 10.^r; %valeur des nouveaux points en normal
        %    sup = 10.^sup; %valeur maximale des nouveaux intervalles en normal
        %    epsi = sup - value; %epsi des nouveaux intervalles en normal
        
        %possibility 3 : define value using the lower and higher bounds
        %    inf = r - epsi; %valeur minimal des nouveaux intervalles en log
        %    sup = r + epsi; %valeur maximale des nouveaux intervalles en log
        %    inf = 10.^inf; %valeur minimal des nouveaux intervalles en normal
        %    sup = 10.^sup; %valeur maximale des nouveaux intervalles en normal
        %    epsi = (sup-inf)/2; %epsi des nouveaux intervalles en normal
        %    value = inf + epsi; %valeur des nouveaux points en normal
        
        position = FindParam(P,params{i}); %on cherche le No de ligne du paramètre
        PRLog.pts(position,indicePts:indicePts+N-1) = value;
        position = PRLog.dim==position; % on cherche le No de ligne de l'epsi
        PRLog.epsi(position,indicePts:indicePts+N-1) = epsi;
    end
    
    indicePts = indicePts + N;
end


if (isfield(P,'selected'))
    PRLog.selected = zeros(1, size(PRLog.pts,2));
end

%TODO : régler les histoires de traj

end
