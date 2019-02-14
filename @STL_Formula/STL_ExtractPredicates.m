function [mus, signals] = STL_ExtractPredicates(phi)
%STL_EXTRACTPREDICATES extracts all predicates and signals present in phi
%
% Synopsis: mus = STL_ExtractPredicates(phi)
%
% Input:
%  - phi : the formula from which the predicates are extraced
%
% Output:
%  - mus: the list of predicates in phi. If phi is a predicate, mus equals
%          phi
%  - signals: list of signals involved in phi (note: currently only detec-
%    ted from patterns of the form signal_id[t])
%
mus = [];
ids = {}; 
if strcmp(phi.type, 'predicate')
    mus = phi;
else
    if ~isempty(phi.phi)
        if ~ismember(phi.id, ids)
            mus = [mus, STL_ExtractPredicates(phi.phi)];
        end
    end
    if ~isempty(phi.phi1)
        mus = [mus, STL_ExtractPredicates(phi.phi1)];
    end
    if ~isempty(phi.phi2)
        mus = [mus, STL_ExtractPredicates(phi.phi2)];
    end
    
    if ~isempty(phi.phin)
        for ii=1:numel(phi.phin)
            mus = [mus, STL_ExtractPredicates(phi.phin(ii))]; %#ok<AGROW>
        end
    end
end

ids = {};
for  ni= 1:numel(mus)
    ids{ni} = mus(ni).id;
end
[~, idx] = unique(ids);
mus = mus(idx);

signals = {};
for imu=1:numel(mus)
     fn_ = mus(imu).params.fn;
     [~,~, ~, matches, tokens] = regexp(fn_, '(\<\w+\>)[.+?\]');
     for im=1:numel(matches)
        signals{end+1} = tokens{im}{1};
     end          
end
if ~isempty(signals)
   signals = unique(signals); 
end

end
