%% 1. Creating log "allLogFalsify"
display('******** Start tuto_falsify_w_disk_caching.m ********')
tuto_falsify_w_disk_caching;

%% 2. Use-case is storing log on different file server 
display('******** Rename CachingFolder Name ********')
movefile('allLogFalsify','fileServerFolder/allLogFalsify')
clear;

%% 3. Use stored folder again
display('******** Recreate Blog ********')
falsif_pb = FalsificationProblem.load_runs('fileServerFolder/allLogFalsify');
Blog = falsif_pb.GetLog();
BreachSamplesPlot(Blog);