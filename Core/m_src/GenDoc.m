function [success, msg, msg_id, log, index_path] = GenDoc(list_scripts, varargin)
%  GenDoc takes a list of scripts and create html documentation

InitBreach;

% Default list of scripts
if nargin==0||isempty(list_scripts)
    list_scripts = get_BrDemo_scripts();
end

% Default publish dir
publish_dir = [get_breach_path( ) filesep 'Doc' ];

options = struct('format', 'html', 'PublishDir', publish_dir, 'log_file', []);
if ~isempty(varargin)
    options = varargin2struct_breach(options, varargin{:});
end

docrun_map = CreateDocRunMap(list_scripts, options.PublishDir);

% Not sure we need the following
if ~isempty(options.log_file)
    [fid, msg] = fopen( options.log_file, 'w');
    if isequal(fid,-1)
        error(['Couldn''t open log file' log_file '. ' msg]);
    end
else
    fid =[];
end

if ischar(list_scripts)
    list_scripts = { list_script };
end

%% if html, create index
if isequal(options.format, 'html')
    index_path = CreateHTMLIndex(list_scripts,  options.PublishDir);
end

[success, msg, msg_id, log] = PublishAll(docrun_map, list_scripts, fid, options.format);
close all;
% create failed log if necessary
i_failed = find(success~=true);
if ~isempty(i_failed)&&~isempty(fid)
    log_failed =  ['FAILED_' log_file];
    fid_failed = fopen(log_failed,'w');
    for i= i_failed
        log_msg = [list_scripts{i} ' FAILED with error msg id:' msg_id{i} '\n'];
        fprintf(fid_failed, log_msg );
    end
    fclose(fid_failed);
end

end

function DocRunMap = CreateDocRunMap(list_scripts, publish_dir)
DocRunMap = containers.Map();
num_scripts = numel(list_scripts);
for i_sc = 1:num_scripts
    name= list_scripts{i_sc};
    dc = DocRun(name);
    % Create dir
    dc.publish_dir = publish_dir;
    dc.publish_html_dir = [publish_dir filesep 'html'];
    
    [success, msg] = mkdir(dc.publish_dir);
    if success==0
        error(msg);
    end
    [success, msg] = mkdir(dc.publish_html_dir);
    if success==0
        error(msg);
    end
    
    DocRunMap(name) = dc;
    
end

end

function [success, msg, msg_id, log] = PublishAll(docrun_map, list_scripts, fid, format)
%

num_flds = numel(list_scripts);
log = '\n--------------------------------------------------------------------\n';
fprintf(log);

for i_sc = 1:num_flds
    sc = list_scripts{i_sc};
    try
        dc = docrun_map(sc);
        switch format
            case 'beamer'
                dc.publish_beamer(1,0);
            case 'html'
                dc.publish_html();
            case 'none'
                dc.run();
        end
        
        success(i_sc) = true;
        msg{i_sc} = '';
        msg_id{i_sc} = '';
    catch
        success(i_sc)= false;
        [msg{i_sc}, msg_id{i_sc}] = lasterr;
    end
    if success(i_sc)
        log_msg = ['***PASSED***  ' regexprep(list_scripts{i_sc},'\\', '\\\')  '\n'];
    else
        log_msg = ['***FAILED****  ' regexprep(list_scripts{i_sc},'\\', '\\\') '. Error msg:' msg{i_sc} '\n'];
    end
    if ~isempty(fid)
        fprintf(fid, log_msg);
    end
    log = [ log log_msg];
    fprintf(log_msg);
end
log = [log '--------------------------------------------------------------------\n'];

fprintf('--------------------------------------------------------------------\n');
if ~isempty(fid)
    fclose(fid);
end

end

function index_path = CreateHTMLIndex(list_script, publish_dir)
% CreateHTMLIndex

index_tmp_file = 'index.m';
fid = fopen( index_tmp_file , 'w+');
if isequal(fid,-1)
    error(['Couldn''t create index file. ' msg]);
end

fprintf(fid, '%%%% List of scripts\n%%\n');
fprintf(fid, '%% Generated on %s with Breach Version %s\n', char(datetime), BreachVersion());
fprintf(fid, '%%%%%%\n');
fprintf(fid, '%% <html>\n%% <ul>\n');
for is =1:numel(list_script)
    [sc_path, sc_name]  = fileparts(which(list_script{is}));
    web_ref = ['%% <li> <a href=" html/'  ...
        sc_name '.html">'...
        sc_name '</a></li>\n'];
    fprintf(fid, web_ref);
end
fprintf(fid, '%% </ul>\n%% </html>');
fclose(fid);
index_path = publish(index_tmp_file, 'outputDir', publish_dir);
delete(index_tmp_file);
disp([ 'Index file created at <a href="matlab: web(''' index_path ''')"> ' index_path '</a>'] );

end

function list_brdemo_scripts = get_BrDemo_scripts()
% TODO discover stuff

BrDemo_path = [get_breach_path() filesep 'Examples' filesep 'BrDemo'];
% list scripts in there 

% Autotrans tuto 
list_brdemo_scripts = {'BrDemo.Autotrans_tutorial'};

% Abstract Fuel Control model
list_AFC_demos = {...
    'BrDemo.AFC_1_Interface',...
    'BrDemo.AFC_2_Simulation',...
    'BrDemo.AFC_3_Sets',...
    'BrDemo.AFC_4_Specifications',...
    'BrDemo.AFC_5_1_Falsification',...
    'BrDemo.AFC_5_2_ParamSynth',...
    'BrDemo.AFC_5_3_MaxSat',...
    'BrDemo.AFC_5_4_ReqMining',...
    };

% AFC Online
list_AFC_demos = [list_AFC_demos 'BrDemo.AFC_Online_Monitoring'];

list_brdemo_scripts = [ list_AFC_demos list_brdemo_scripts ];



end


