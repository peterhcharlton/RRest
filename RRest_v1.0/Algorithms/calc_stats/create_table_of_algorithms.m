function create_table_of_algorithms(up)
%CREATE_TABLE_OF_ALGORITHMS creates a table of algorithms and their performances
%	            create_table_of_algorithms(up)
%
%	Inputs:
%		data            data files should be stored in the specified format
%       up              universal parameters structure
%
%	Outputs:
%       tables          csv tables of results
%

fprintf('\n--- Creating Results Tables ');

% save setup
save_name = up.paths.filenames.results_table;
savepath = [up.paths.tables_save_folder, save_name, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

% Load BA Data
load_name = 'BA_results';
loadpath = [up.paths.data_save_folder, up.paths.filenames.global_BA, '.mat'];
load(loadpath, load_name);

% Load algorithm names
load_name = up.paths.filenames.alg_names;
loadpath = [up.paths.data_save_folder, up.paths.filenames.alg_names, '.mat'];
load(loadpath, load_name);

% Load study data
load_name = 'trad_stats';
loadpath = [up.paths.data_save_folder, up.paths.filenames.study_stats, '.mat'];
load(loadpath, load_name);

% identify groups
groups = fieldnames(BA_results);

%% Create a BA table and a Stats table for each group

for group_no = 1 : length(groups)
    group = groups{group_no};
    
    %% Extract relevant data
    rel_data = extract_rel_data(group, BA_results, trad_stats, up);
    
    %% Make Table of BA Results
    BA_results_table = create_table_results(alg_names, rel_data, up, group, 'BA');
    
    %% Save to file
    save_name = [up.paths.filenames.BA_results_table, '_', group];
    savepath = [up.paths.tables_save_folder, save_name];
    xlswrite(savepath, BA_results_table)
    
    %% Make Table of Stats Results
    stats_results_table = create_table_results(alg_names, rel_data, up, group, 'stats');
    
    %% Save to file
    save_name = [up.paths.filenames.stats_results_table, '_', group];
    savepath = [up.paths.tables_save_folder, save_name];
    xlswrite(savepath, stats_results_table)
    
end

end

function rel_data = extract_rel_data(group, BA_results, stats, up)

rel_data.BA = extractfield(BA_results, group); rel_data.BA = rel_data.BA{1,1};
rel_data.stats = extractfield(stats, group); rel_data.stats = rel_data.stats{1,1};

end

function table_data = create_table_results(alg_names, rel_data, up, group, type)

% Find alg numbers:
alg_names.no = cell(length(alg_names.names),1);
for s = 1 : length(alg_names.no)
    alg_names.no{s} = ['a', num2str(s)];
end

% Extract relevant variables for table:
tab_vars.alg_no = alg_names.no;
if isfield(alg_names.meths, 'xa')
    tab_vars.m_xa = alg_names.meths.xa;
    tab_vars.m_xb = alg_names.meths.xb;
    tab_vars.m_ef = alg_names.meths.ef;
    tab_vars.m_et = alg_names.meths.et;
    tab_vars.m_fm = alg_names.meths.fm;
    tab_vars.m_ft = alg_names.meths.ft;
    % modify the fme:
    rel_els = ~isnan(tab_vars.m_fm);
    tab_vars.m_xb(rel_els) = 99;
end
tab_vars.alg_name = alg_names.names;
tab_vars.signal = alg_names.sigs;

% Extract relevant statistics for table:
if strcmp(type, 'BA')
    tab_vars.bias = rel_data.BA.bias.val;
    tab_vars.two_sd = rel_data.BA.two_sd.val;
    tab_vars.prop_wins_est = rel_data.stats.prop_wins_est;
    tab_vars.prop_wins_ref_and_good_sqi = rel_data.stats.prop_wins_good_sqi_and_ref;
elseif strcmp(type, 'stats')
    tab_vars.percerr = rel_data.stats.percerr;
    tab_vars.mae = rel_data.stats.mae;
    tab_vars.sdae = rel_data.stats.sdae;
    tab_vars.rmse = rel_data.stats.rmse;
    tab_vars.prop_wins_est = rel_data.stats.prop_wins_est;
end

%% Determine what goes into the table
headers = fieldnames(tab_vars);
table_data = cell(1+length(tab_vars.alg_name), length(headers));
table_data(1,:) = headers;
for col_no = 1 : length(headers)
    eval(['col_data = tab_vars.' headers{col_no} ';']);
    if ~iscell(col_data)
        table_data(2:end,col_no) = num2cell(col_data);
    else
        table_data(2:end,col_no) = col_data;
    end
end

end


function [subj_data, BA_data, stats_data, old_names] = load_data(group, up)

% Load subj data
load_name = ['data_', group];
loadpath = [up.paths.data_save_folder, up.paths.filenames.group_data, group, '.mat'];
load(loadpath, load_name)
eval(['subj_data = data_' group  ';']);

% Load B-A results
load_name = 'BA';
loadpath = [up.paths.data_save_folder, up.paths.filenames.group_BA, group, '.mat'];
load(loadpath, load_name)
BA_data = BA; clear BA

% Load Stats data
load_name = 'acc';
loadpath = [up.paths.data_save_folder, up.paths.filenames.group_stats, group, '.mat'];
load(loadpath, load_name)
stats_data = acc;
old_names = acc.rr_alg_names;

end

function create_table_stats_results_subj(alg_names, sig_names, old_names, subj_data, up, group)

%% Find prop with an abs error of < 1 bpm
eval(['subjs = unique(subj_data.' old_names{1} '.subj);']);
prop_acc = zeros(length(old_names), length(subjs));
for alg_no = 1 : length(old_names)
     eval(['rel_data = subj_data.' old_names{alg_no} ';']);
     
     rel_data.error = rel_data.est - rel_data.ref;
     rel_data.abs_error = abs(rel_data.error);
     rel_data.acc_log = rel_data.abs_error < 1;
     
     for subj_no = 1:length(subjs)
         subj = subjs(subj_no);
         rel_els = logical(rel_data.subj == subj);
         prop_acc(alg_no, subj_no) = sum(rel_data.acc_log(rel_els)) / length(rel_data.acc_log(rel_els));
     end
     
     % extract relevant parts of algs:
     for comp_str = {'FMe', 'Flt', 'RSa', 'PDt', 'Est', 'MFu', 'TFu'}
         if ~exist([lower(comp_str{1,1}), '_meth'], 'var')
             eval([lower(comp_str{1,1}), '_meth = nan(length(old_names),1);']);
         end
         rel_el = strfind(alg_names{alg_no}, comp_str{1,1});
         if rel_el
             temp = str2num(alg_names{alg_no}(rel_el+3));
             eval([lower(comp_str{1,1}), '_meth(alg_no) = temp;']);
         end
     end
     
     clear rel_data
end


%% Determine what goes into the table
headers = {'Signal', 'Algorithm Name', 'Old Name', 'FMe', 'Flt', 'RSa', 'PDt', 'Est', 'MFu', 'TFu'};
for subj_no = 1:length(subjs)
    headers{end+1} = ['Prop_acc ' num2str(subjs(subj_no))];
end


table_data = cell(1+length(alg_names), length(headers));
table_data(1,:) = headers;
table_data(2:end,1) = sig_names;
table_data(2:end,2) = alg_names;
table_data(2:end,3) = old_names;
table_data(2:end,4) = num2cell(fme_meth);
table_data(2:end,5) = num2cell(flt_meth);
table_data(2:end,6) = num2cell(rsa_meth);
table_data(2:end,7) = num2cell(pdt_meth);
table_data(2:end,8) = num2cell(est_meth);
table_data(2:end,9) = num2cell(mfu_meth);
table_data(2:end,10) = num2cell(tfu_meth);
for s = 1 : size(prop_acc,2)
    table_data(2:end,10+s) = num2cell(prop_acc(:,s));
end

%% Save to file
save_name = [up.paths.filenames.subj_results_table '_' group];
savepath = [up.paths.tables_save_folder, save_name];
xlswrite(savepath,table_data)

end

function create_table_BA_results(alg_names, sig_names, BA_data, up, group)

%% Determine what goes into the table
headers = {'Signal', 'Algorithm Name', 'Bias', 'Bias (LCI)', 'Bias (UCI)', 'LLOA', 'LLOA (LCI)', 'LLOA (UCI)', 'ULOA', 'ULOA (LCI)', 'ULOA (UCI)', 'Precision', 'Precision (LCI)', 'Precision (UCI)'};
table_data = cell(1+length(alg_names), length(headers));
table_data(1,:) = headers;
table_data(2:end,1) = sig_names;
table_data(2:end,2) = alg_names;
table_data(2:end,3) = num2cell(BA_data.bias.val);
table_data(2:end,4) = num2cell(BA_data.bias.lci);
table_data(2:end,5) = num2cell(BA_data.bias.uci);
table_data(2:end,6) = num2cell(BA_data.lloa.val);
table_data(2:end,7) = num2cell(BA_data.lloa.lci);
table_data(2:end,8) = num2cell(BA_data.lloa.uci);
table_data(2:end,9) = num2cell(BA_data.uloa.val);
table_data(2:end,10) = num2cell(BA_data.uloa.lci);
table_data(2:end,11) = num2cell(BA_data.uloa.uci);
table_data(2:end,12) = num2cell(BA_data.prec.val);
table_data(2:end,13) = num2cell(BA_data.prec.lci);
table_data(2:end,14) = num2cell(BA_data.prec.uci);

%% Save to file
save_name = [up.paths.filenames.BA_results_table '_' group];
savepath = [up.paths.tables_save_folder, save_name];
xlswrite(savepath,table_data)

end

function create_table_stats_results(alg_names, sig_names, stats_data, up, group)

%% Determine what goes into the table
headers = {'Signal', 'Algorithm Name', 'MAE', 'MAE (lq)', 'MAE (uq)', 'SDAE', 'SDAE (lq)', 'SDAE (uq)', 'PERC ERR', 'PERC ERR (lq)', 'PERC ERR (uq)', 'RMSE', 'RMSE (lq)', 'RMSE (uq)', 'MEAN ERROR', 'MEAN ERROR (lq)', 'MEAN ERROR (uq)', 'STD ERROR', 'STD ERROR (lq)', 'STD ERROR (uq)', 'PROP WINS', 'PROP WINS (lq)', 'PROP WINS (uq)'};
table_data = cell(1+length(alg_names), length(headers));
table_data(1,:) = headers;
table_data(2:end,1) = sig_names;
table_data(2:end,2) = alg_names;
table_data(2:end,3) = num2cell(stats_data.mae.med);
table_data(2:end,4) = num2cell(stats_data.mae.lq);
table_data(2:end,5) = num2cell(stats_data.mae.uq);
table_data(2:end,6) = num2cell(stats_data.sdae.med);
table_data(2:end,7) = num2cell(stats_data.sdae.lq);
table_data(2:end,8) = num2cell(stats_data.sdae.uq);
table_data(2:end,9) = num2cell(stats_data.percerr.med);
table_data(2:end,10) = num2cell(stats_data.percerr.lq);
table_data(2:end,11) = num2cell(stats_data.percerr.uq);
table_data(2:end,12) = num2cell(stats_data.rmse.med);
table_data(2:end,13) = num2cell(stats_data.rmse.lq);
table_data(2:end,14) = num2cell(stats_data.rmse.uq);
table_data(2:end,15) = num2cell(stats_data.mean_diff.med);
table_data(2:end,16) = num2cell(stats_data.mean_diff.lq);
table_data(2:end,17) = num2cell(stats_data.mean_diff.uq);
table_data(2:end,18) = num2cell(stats_data.std_diff.med);
table_data(2:end,19) = num2cell(stats_data.std_diff.lq);
table_data(2:end,20) = num2cell(stats_data.std_diff.uq);
table_data(2:end,21) = num2cell(stats_data.prop_wins.med);
table_data(2:end,22) = num2cell(stats_data.prop_wins.lq);
table_data(2:end,23) = num2cell(stats_data.prop_wins.uq);

%% Save to file
save_name = [up.paths.filenames.stats_results_table '_' group];
savepath = [up.paths.tables_save_folder, save_name];
xlswrite(savepath,table_data)

end