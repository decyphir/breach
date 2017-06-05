function CleanModelsData
% Delete files created by Breach, model copies, etc

mdl_data_path = [ get_breach_path() filesep 'Ext' filesep 'ModelsData' ];
ls(mdl_data_path)
go = input('Delete all files?');
if go 
   delete([mdl_data_path filesep '*']);
end
 