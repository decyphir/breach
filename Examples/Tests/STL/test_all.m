function Reports = test_all()
% test_all run and check all scripts
% TODO check status for passed or failed status
% TODO more generally improve output
% TODO add more test scripts

%% Init Breach
clear
InitBreach;
Reports = {};

%% Test_script1
T1 = DocRun('test_script1');
R1 = T1.run_and_cmp();
printAll(R1);
Reports = [Reports R1]; 

%% Add more test scripts..


end

function printAll(R)
    fn = fieldnames(R);
    
    for i_f = 1:numel(R)
        fprintf(['Result for variable ' fn{i_f} ':\n'])
        R.(fn{i_f}).printStatus();
        fprintf('\n');
    end
end