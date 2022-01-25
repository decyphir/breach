function Nn = N2Nn(n,nb, num_new, varargin)
%N2NN Builds a set of nb^n points in N^n (nb points per axes)
%
% Synopsis : Nn = N2Nn(n,nb, [num_new, varargin])
%
%  The resulting set of points is sorted according to a heuristic which favors clusters of identical values.
%  In particular, the two first column are uniform. In the corner
%  computation case, this corresponds to the two extremal corners.
%

opt.randomize_order = true;
opt.start_with_min_max = true;
opt.sort_by_clusters = false;
opt = varargin2struct_breach(opt,varargin{:});

if isscalar(nb)
    nb = nb*ones(1,n);
end

global BreachGlobOpt

max_num_new = prod(nb);
if nargin<=2
    num_new = max_num_new;
else
    num_new = min(max_num_new, num_new);
end

if num_new >BreachGlobOpt.MaxNumSamples
    warning('Maximum number %d of samples reached. Change global variable BreachGlobOpt.MaxNumSamples if you need more.',    BreachGlobOpt.MaxNumSamples);
    num_new =BreachGlobOpt.MaxNumSamples;
end

if n==1
    Nn = 1:num_new;
elseif num_new ==1
    Nn= ones(n,1);
elseif num_new ==2
    Nn=  [ ones(n,1) nb'];
else
    p = 1:nb(1);
    if num_new == max_num_new
        Nn = N2NnIter(n,p,nb);
    else % this is where we need some smarts, let's do simple first
        rng(1, 'twister');  % seed random generator to get deterministic results
        if max_num_new < 1e5 %  if max_num_new is still reasonable, just truncate full grid
            Nn =  N2NnIter(n,p,nb);
            %FullNn =  N2NnIter(n,p,nb);
            %idx= randperm(max_num_new);
            %Nn = FullNn(:, idx(1:num_new+1));
            % always include/starts with min and max corners
            %Nn =  unique([ones(n,1) nb' Nn(:,2:end)]', 'rows','stable')';
        else % otherwise random stuff until we get enough unique vectors
            Nn = unique(RandNn(n,nb,1e6)', 'rows', 'stable')';
            while size(Nn, 2) < num_new
                Nn_more = RandNn(n,nb,num_new);
                Nn =  unique([Nn Nn_more]', 'rows','stable')';
            end
            % always include/starts with min and max corners
            Nn =  unique([ones(n,1) nb' Nn(:,2:end)]', 'rows', 'stable')';
        end
    end
    % optimize order
    % max min ( size(biggest_cluster of zeros, biggest_cluster_of_ones), )
    clust = find_size_min_cluster(Nn);
    [~, clust_sort ] = sort(clust, 1, 'descend');
    Nn = Nn(:, clust_sort);
    
    % throw in max frequency in 3 and 4 for good measure
    h = size(Nn,1);
    onetwos = repmat([1 2]',ceil(h/2),1);
    onetwos = onetwos(1:h,1);
    twoones = repmat([2 1]',ceil(h/2),1);
    twoones = twoones(1:size(Nn,1),1);       
    Nn = unique([Nn(:,1:2) twoones onetwos Nn(:,3:end)]', 'rows', 'stable')';        
    Nn = Nn(:,1:num_new);
end
end

function num = find_size_min_cluster(x)
%
x = diff(x)';
num  = zeros(size(x,1),1);
for irow = 1:size(x, 1)
    i_non_zero =[0 find(x(irow,:)) size(x,2)+1];
    num_non_zeros =  diff(i_non_zero)-1;
    num(irow,1)= min(num_non_zeros);
end

end

function notNn = RandNn(n, nb, num_new)
% Generate random in the grid
notNn = zeros(n, num_new);
for in = 1:n
    notNn(in,:) = randi(nb(in),1, num_new);
end

end

function Nn = N2NnIter(n,p, nb)
if(n==1)
    Nn = p;
else
    Nnm = N2NnIter(n-1,p,nb(1:end-1));
    nbm = size(Nnm,2);    
    Nn = [ repmat(Nnm,1,nb(end)) ; kron(1:nb(end),ones(1,nbm)) ];
end
end
