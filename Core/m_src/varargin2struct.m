function strct  = varargin2struct(strct, varargin)
% varargin2struct(default_strct, varargin) Convert varargin given as pair
% of name, values into a struct. Default is provided as strct.

% read the acceptable names
optionNames = fieldnames(strct);

% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('Options should be provided as  propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   inpName = pair{1}; %# make case insensitive
   idx =strcmpi(inpName,optionNames);
   if any(idx)
      %# overwrite options. If you want you can test for the right class here
      %# Also, if you find out that there is an option you keep getting wrong,
      %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
      strct.(optionNames{idx}) = pair{2};
   else
      error('%s is not a recognized option name',pair{1})
   end
end

end
