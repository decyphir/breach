function U = UniStepSimulinkInput(cp, Sys, pts, tspan)
% currently obsolete - to be removed soon 


if isstruct(Sys)
    DimU = Sys.DimU;
else
    InputNames = Sys; % not very nice.
end


if (isempty(pts)&&isempty(tspan))
    
    U.params = {};
    U.p0 = [];
    
    for ku = 1:numel(InputNames)
        for k = 1:cp
            U.params = {U.params{:} [InputNames{ku} '_u' num2str(k-1)]};
            U.p0 = [U.p0 0];
        end
    end
else
    
    U.t = 0;
    U.u = reshape(pts(end-cp*DimU+1:end), [cp DimU] );
    U.u = [U.u ; U.u(end,:)];
    
    for k=2:cp
        U.t(k,1) = (k-1)*tspan(end)/cp;
    end
    U.t(cp+1,1) = tspan(end);
    
    U.u = reshape( [U.u'; U.u'], [size(U.u',1) 2*size(U.u',2)] )';
    U.t = reshape( [U.t'; [U.t(2:end)' U.t(end)'*2]-10*eps], [1 2*numel(U.t)] )';
    
end
