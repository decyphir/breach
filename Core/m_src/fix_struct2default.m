function s  = fix_struct2default(default, s)
% fix_struct2defaut(default, s) copy default field if absent in s

fs = fieldnames(default);
for ii= 1:numel(fs)
   if ~isfield(s,fs{ii}) 
      s.(fs{ii})=default.(fs{ii});
   end   
end
