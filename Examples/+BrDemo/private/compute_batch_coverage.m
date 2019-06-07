function [covocc, covlog, covent] = compute_batch_coverage(B, num_samples_array,  epsi, delta)

covlog = [];
covocc = []; 
covent = [];
num_var = numel(B.GetVariables());
for i= num_samples_array
    Bi = B.copy();
    Bi.SetEpsGridsize(epsi*ones(1,num_var));
    Bi.SetDeltaGridsize(delta*ones(1,num_var));
    Bi.QuasiRandomSample(i);
    Bi.AddPoints();
    
    covocc= [covocc Bi.ComputeCellOccupancyCoverage()];
    covlog = [covlog Bi.ComputeLogCellOccupancyCoverage()]; 
    covent = [covent Bi.ComputeEntropyCoverage()]; 
    
end

end