function st_err = ME2st(ME)
st_err = sprintf('Error id: %s\nError message: %s', ME.identifier, ME.message);
rc = sprintf('\n');
for i = 1:size(ME.stack)
    st_err = [st_err rc 'file: ' ME.stack(i,1).file rc ME.stack(i,1).name ' line:' num2str(ME.stack(i,1).line)];
end
end
