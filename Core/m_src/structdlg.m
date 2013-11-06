function [P, answer] = structdlg(P,title,prompt, default)
  
 % helper function that opens a dialog window to edit the fields of a
 % given structure
   
   if (nargin == 1)
     title = 'Edit params';
   end
   if (~exist('prompt'))
     prompt= fieldnames(P);
   end
   
   if (~exist('default'))
     default_values = struct2cell(P);
   elseif (isempty(default))
       default_values = struct2cell(P);
   else
       default_values = default;
   end     
     
   for i = 1:numel(prompt)
     if isnumeric(default_values{i})
       defaultanswer{i}= num2str(default_values{i});
     elseif ~isstr(default_values{i}) 
       defaultanswer{i}= '';
     else
       defaultanswer{i}=default_values{i};
     end
   end
         
   options.Resize='on';
   options.WindowStyle='normal';
   options.Interpreter='tex';
 
   answer = inputdlg(prompt,title ,1,defaultanswer,options); 
   if isempty(answer)
       P = [];
       return
   end
   fnames = fieldnames(P);
   for i = 1:numel(prompt)
       %if (isstr(default_values{i}))
       %  P.(fnames{i}) = answer{i}; 
       %  continue;
       %end
       if (~strcmp(answer{i},'')&&~(isempty(answer{i})))
         try 
           P.(fnames{i}) = eval(answer{i});
         catch
           warndlg('Invalid value');
         end
       end
   end