function U = VarStepSimulinkInput(cp, Sys, pts, tspan)  
% FIXME this function does two things, which is bad: it names input signals or compute input signals
% It is currently obsolete - to be removed soon 
  
if isstruct(Sys)
    DimU = Sys.DimU;
else
    InputNames = Sys; % not very nice.
end

  if (isempty(pts)&&isempty(tspan))

    U.params = {};
    U.p0 = [];

    for k = 1:cp    
      U.params = {U.params{:}  ['dt_u' num2str(k-1)]};
      U.p0 = [U.p0 1];
      for ku = 1:numel(InputNames)
        U.params = {U.params{:} [InputNames{ku} '_u' num2str(k-1)]};
        U.p0 = [U.p0 0];
      end    
    end
  
  else    
                
    U.u = reshape(pts(end-cp*(DimU+1)+1:end), [(DimU+1) cp] )';
    U.t =[0 ; cumsum(U.u(:,1))];
    U.u = U.u(:,2:end);
    U.u = [U.u ; U.u(end,:)];        

  end
