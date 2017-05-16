function [B, summary, traces] = LoadResults(folder_name, verbose)
% LoadResults(folder_name) Import results saved previously with SaveResults command 
%
% FIXME Overall a lazy implementation, a bit too much of redundancy

if ~exist('verbose', 'var')
    verbose = 1;
end

%% checks whether there is a breach system 
breach_system_fn = [folder_name filesep 'breach_system.mat' ];
if ~exist(breach_system_fn, 'file')
    error('Breach:LoadResult:FileNotFound', 'Result files not found.');
end

ss = load(breach_system_fn); 
Bname = fieldnames(ss);
Bname = Bname{1};
B = ss.(Bname);
[summary, traces] = B.ExportTracesToStruct();

switch nargout
    case 0
      assignin('base', Bname, ss.(Bname));
      if verbose
          disp(['Loaded ' Bname  ' from ' breach_system_fn '.']);
      end
    case 1
       B = ss.(Bname);       
end