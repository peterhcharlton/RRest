function calc_stats(up)
%CALC_STATS calculates statistics based on the estimated and ref RRs
%	            calc_stats(up)
%
%	Inputs:
%		data            data files should be stored in the specified format
%       up              universal parameters structure
%
%	Outputs:
%       win_data        files containing a summary of each window analysed
%       BA              files containing Bland-Altman analysis results
%       study_stats     file containing other statistical analysis results
%       alg_names       file containing algorithm names
%

%% Calculate stats for impedance numerics if in dataset
if up.analysis.imp_stats
    calc_stats_for_impedance(up);
end

%% Decide whether this is a single period, or a merged period

if isempty(strfind(lower(up.paths.data_load_filename), 'rest_and_rec'))
    % then it is a single period
    calc_stats_for_a_single_period(up);
else
    % then it is a merged period
    calc_stats_for_rest_and_rec(up);
end

if ~isempty(strfind(lower(up.paths.data_load_filename), 'vortal_rest_and_rec_data_clin'))
    create_table_win_data_merged_raw_clin_imp(up);
    create_table_win_data_merged_raw_clin(up);
end


end

%% Impedance scripts

function calc_stats_for_impedance(up)

if isempty(strfind(lower(up.paths.data_load_filename), 'rest_and_rec'))
    %% Create table of each algorithm, each subject, and each window
    create_table_win_data_imp(up);
    
    %% Calculate stats for entire study, and for each sub-group analysis
    calc_study_and_sub_group_stats_imp(up);
else
    % then it is a merged period
    %% Create table of each algorithm, each subject, and each window
    create_table_win_data_merged_imp(up);
    
    %% Calculate stats for entire study, and for each sub-group analysis
    calc_study_and_sub_group_stats_imp(up);
end

end

function create_table_win_data_imp(up)

%% save setup
save_name = up.paths.filenames.win_data;
savepath = [up.paths.data_save_folder, up.paths.filenames.win_data_imp, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

fprintf('\n--- Creating a Table of Impedance Results for each window ');

%% Specify algorithm names
alg_names.old_names = {'num'};
alg_names.names = {'Imp_num'};

%% Find out which group each subject is in
subjs = up.paramSet.subj_list;
rel_els = find_rel_subjs(up);
groups = fieldnames(rel_els);
for group_no = 1 : length(groups)
    eval([groups{group_no} '_log = false(length(subjs),1);']);
    eval([groups{group_no} '_log(rel_els.' groups{group_no} ') = true;']);
end

%% Cycle through each algorithm, subject and window
ecg_log = true(length(alg_names.old_names),1);        % is it ecg or ppg

% create storage for raw data
[data.subj, data.ecg_log, data.alg_no] = deal(single([]));
for group_no = 1 : length(groups)
    eval(['data.' groups{group_no} '_log = single([]);']);
end
[data.win_no, data.est, data.ref, data.sqi, data.hr, data.hr_rr, data.snr_log] = deal(single([]));
data.alg_names = cell(0);

% load ref RRs and est RRs for all algorithms for all subjects
counter_no = 0;
for subj = subjs
    
    % Load ref RRs
    loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.ref_rrs, '.mat'];
    load(loadpath);
    
    % Load est RRs from imp
    loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.ref_rrs, '.mat'];
    temp = load(loadpath);
    est_RRs = rr_ref.imp;
    
    % Load SQIs
    loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.sqi, '.mat'];
    load(loadpath);
    
    % Store data
    ppg_hr_rr_ratio = sqi.ppg.hr./rr_ref.v;
    ekg_hr_rr_ratio = sqi.ekg.hr./rr_ref.v;
    for s = 1 : length(alg_names.old_names)
        rel_data = est_RRs.v(:);
        start_el = counter_no+1;
        counter_no = counter_no + length(rel_data);
        data.est(start_el:counter_no,1) = rel_data;
        data.alg_no(start_el:counter_no,1) = s;
        if ecg_log(s)
            data.sqi(start_el:counter_no,1) = sqi.ekg.v;
            data.hr(start_el:counter_no,1) = sqi.ekg.hr;
            data.hr_rr(start_el:counter_no,1) = ekg_hr_rr_ratio;
        else
            data.sqi(start_el:counter_no,1) = sqi.ppg.v;
            data.hr(start_el:counter_no,1) = sqi.ppg.hr;
            data.hr_rr(start_el:counter_no,1) = ppg_hr_rr_ratio;
        end
    end
    
    % variables which are already in the order of the windows
    data.ref = [data.ref; repmat(rr_ref.v(:), [length(alg_names.old_names),1])];
    data.snr_log = [data.snr_log; repmat(rr_ref.snr_log(:), [length(alg_names.old_names),1])];
    win_nos = 1 : length(rr_ref.v);
    data.win_no = [data.win_no; repmat(win_nos(:), [length(alg_names.old_names),1])];
    
    % variables which are in the order of the algorithms , and need to be changed
    temp = repmat(alg_names.names(:)', [length(rr_ref.v), 1]);
    temp = temp(1:numel(temp)); temp = temp(:);
    data.alg_names = [data.alg_names; temp];
    temp = repmat(ecg_log(:)', [length(rr_ref.v), 1]);
    temp = temp(1:numel(temp)); temp = temp(:);
    data.ecg_log = [data.ecg_log; temp];
    
    % variables which are fixed for this subj
    data.subj = [data.subj; subj*ones(length(rr_ref.v)*length(alg_names.old_names),1)];
    %data.young_log = [data.young_log; young_log(subj)*ones(length(rr_ref.v)*length(alg_names.old_names),1)];
    for group_no = 1 : length(groups)
        eval(['data.' groups{group_no} '_log = [data.' groups{group_no} '_log; ' groups{group_no} '_log(subj)*ones(length(rr_ref.v)*length(alg_names.old_names),1)];']);
    end
    
end
clear win_nos comp_no s start_el counter_no est_RRs temp rr_ref loadpath subj ecg_log rel_data ppg_hr_rr_ratio sqi ekg_hr_rr_ratio

for group_no = 1 : length(groups)
    eval(['data.' groups{group_no} '_log = logical(data.' groups{group_no} '_log);']);
end
data.ecg_log = logical(data.ecg_log);
data.sqi = logical(data.sqi);
data.snr_log = logical(data.snr_log);
data.comb_log = logical(data.sqi & data.snr_log);

%% Save
algorithm_names = data.alg_names;
data = rmfield(data, 'alg_names');
win_data = data;
save(savepath, save_name, 'algorithm_names')

end

function calc_study_and_sub_group_stats_imp(up)

save_name = 'BA_results';
savepath = [up.paths.data_save_folder, up.paths.filenames.imp_BA, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

fprintf('\n--- Performing BA analyses for entire cohort and sub-groups ');

%% Load data from entire study
load_name = up.paths.filenames.win_data;
loadpath = [up.paths.data_save_folder, up.paths.filenames.win_data_imp, '.mat'];
load(loadpath, load_name)

%% Identify groups of subjects for which to perform BA analyses:
% i) Entire cohort
% ii) Subgroups
% iii) High and Low Quality
% iv) HR (quintiles)
% v) RR (quintiles)
% vi) HR:RR (quintiles)

group_names = unique(up.paramSet.groups);
groups.entire = true(length(win_data.subj),1);
for group_no = 1 : length(group_names)
    eval(['groups.' group_names{group_no} ' = win_data.' group_names{group_no} '_log;']);
end
%groups.low_qual = ~win_data.sqi;

%% remove these groups for the moment
include_extra_groups = 0;

if include_extra_groups
    no_quantiles = 20;
    win_data.rr = win_data.ref;
    for param = {'hr', 'rr', 'hr_rr'}
        eval(['rel_data = win_data.' param{1,1} ';']);
        for quintile_no = 1 : no_quantiles
            lower_boundary = quantile(rel_data, (quintile_no-1)/no_quantiles);
            if quintile_no == 1
                lower_boundary = lower_boundary - 0.1;
            end
            upper_boundary = quantile(rel_data, quintile_no/no_quantiles);
            eval(['groups.' param{1,1} '_' num2str(quintile_no) ' = win_data.' param{1,1} '>lower_boundary & win_data.' param{1,1} ' <= upper_boundary;']);
            clear upper_boundary lower_boundary
        end
        clear rel_data quintile_no
    end
    clear param
end


%% Find errors
win_data.errors = win_data.est - win_data.ref;

%% BA Analyses

% compute for each group
group_names = fieldnames(groups);
% setup
alg_nos = unique(win_data.alg_no);
no_algs = length(alg_nos);

no_subjs = length(unique(win_data.subj));
for group_no = 1 : length(group_names)
    
    group_name = group_names{group_no,1};
    eval(['rel_group_els = groups.' group_name ';']);
    
    [BA.bias.val, BA.bias.lci, BA.bias.uci, ...
        BA.prec.val, BA.prec.lci, BA.prec.uci, ...
        BA.lloa.lci, BA.lloa.val, BA.lloa.uci, ...
        BA.uloa.val, BA.uloa.lci, BA.uloa.uci, ...
        BA.prop] = deal(nan(no_algs,1));
    
    % for each algorithm ...
    
    for alg_no_no = 1 : no_algs
        
        alg_no = alg_nos(alg_no_no);
        
        % identify relevant data
        if strcmp(group_name, 'low_qual')
            rel_els = rel_group_els & win_data.alg_no == alg_no & ~win_data.sqi & win_data.snr_log;
        else
            % otherwise exclude low quality windows
            rel_els = rel_group_els & win_data.alg_no == alg_no & win_data.sqi & win_data.snr_log;
        end
        
        rel_errors = win_data.errors(rel_els);
        rel_subjs = win_data.subj(rel_els);
        subjs = unique(rel_subjs);
        no_subjs = length(subjs);
        
        % calculate prop of wins
        BA.prop(alg_no) = sum(~isnan(rel_errors)/length(rel_errors));
        
        % calculate bias
        BA.bias.val(alg_no) = nanmean(rel_errors);
        
        % calculate bias CI
        n = no_subjs;    % sample size (no subjects)
        doff = n - 1;               % Degrees of Freedom
        tval = tinv(0.975,doff);    % value from T distribution
        std_diffs = std(rel_errors(~isnan(rel_errors)));    % standard dev of diffs
        std_err_diffs = sqrt((std_diffs^2)/n);      % standard err of diffs
        BA.bias.lci(alg_no) = BA.bias.val(alg_no) - (tval*std_err_diffs);
        BA.bias.uci(alg_no) = BA.bias.val(alg_no) + (tval*std_err_diffs);
        
        % calculate loa
        [p,tbl] = anova1(rel_errors,rel_subjs,'off'); % one-way anova of the differences
        MS_subj = tbl{2,4};
        MS_resid = tbl{3,4};
        no_wins = nan(no_subjs,1);
        for subj_no = 1:no_subjs
            no_wins(subj_no) = sum(rel_subjs == subjs(subj_no));
        end
        denominator = ( (sum(no_wins)^2) - sum(no_wins.^2) ) / ( (n-1)*sum(no_wins) );
        heterogeneity_variance = (MS_subj - MS_resid)/denominator;
        total_variance = heterogeneity_variance + MS_resid;
        BA.lloa.val(alg_no) = BA.bias.val(alg_no) - (1.96*sqrt(total_variance));
        BA.uloa.val(alg_no) = BA.bias.val(alg_no) + (1.96*sqrt(total_variance));
        
        % calculate loa CI
        std_err_lim = sqrt((3*(std_diffs^2))/n);      % standard err of the limit of diffs
        BA.lloa.lci(alg_no) = BA.lloa.val(alg_no) - (tval*std_err_lim);
        BA.lloa.uci(alg_no) = BA.lloa.val(alg_no) + (tval*std_err_lim);
        BA.uloa.lci(alg_no) = BA.uloa.val(alg_no) - (tval*std_err_lim);
        BA.uloa.uci(alg_no) = BA.uloa.val(alg_no) + (tval*std_err_lim);
        
        % calculate precision
        BA.prec.val(alg_no) = BA.uloa.val(alg_no) - BA.lloa.val(alg_no);
        BA.prec.lci(alg_no) = BA.prec.val(alg_no) - (2*std_err_lim*tval);
        BA.prec.uci(alg_no) = BA.prec.val(alg_no) + (2*std_err_lim*tval);
        
        % calculate 2SD
        BA.two_sd.val(alg_no) = 2*sqrt(total_variance);
        
    end
    
    eval(['BA_results.' group_name '= BA;']);
end

%% Save
save(savepath, save_name);

end

%% Single period scripts

function calc_stats_for_a_single_period(up)

fprintf('\n--- Calculating Statistics ');

%% Rename algorithms
rename_algorithms(up);

%% Create table of each algorithm, each subject, and each window
create_table_win_data(up);

%% Calculate stats for entire study, and for each sub-group analysis
calc_study_and_sub_group_stats(up);

%% Find specific stats if this is the synthetic dataset
if up.analysis.calc_synth_stats
    specific_synth_stats(up);
end

%% Find specific stats for vortal paper;
calc_vortal_stats = 0;
if calc_vortal_stats
    specific_vortal_stats(up);
    specific_vortal_plots(up);
end

%% Find traditional stats (e.g. MAE):
calc_trad_stats(up);

end

function rename_algorithms(up)

%% save setup
save_name = up.paths.filenames.alg_names;
savepath = [up.paths.data_save_folder, up.paths.filenames.alg_names, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

fprintf('\n--- Renaming Algorithms ');

%% Load algorithm names
loadpath = [up.paths.data_save_folder, '1', up.paths.filenames.rrEsts, '.mat'];
temp = load(loadpath); clear loadpath
old_names = sort(fieldnames(temp));     % used to be: old_names = fieldnames(temp);
alg_names.old_names = old_names;

%% Rename algorithm names according to the nomenclature used in the paper:
[alg_names.names, alg_names.sigs] = deal(cell(length(old_names),1));
for alg_no = 1 : length(old_names)
    
    % setup
    rel_name = old_names{alg_no};
    
    rel_name = strrep(rel_name, 'GCE', 'TPE');
    
    % Eliminate fixed components
    rel_name = strrep(rel_name, 'ELF', '');
    rel_name = strrep(rel_name, 'FPt', '');
    rel_name = strrep(rel_name, 'EHF', '');
    rel_name = strrep(rel_name, 'RDt', '');
    rel_name = strrep(rel_name, 'PDt', '');
    rel_name = strrep(rel_name, 'FMe', '');
    rel_name = strrep(rel_name, 'flt', '');
    
    % Filter-based
    rel_name = strrep(rel_name, 'CCF', 'Xa4');
    rel_name = strrep(rel_name, 'BFi', 'Xa1');
    rel_name = strrep(rel_name, 'Wfm', 'Xa3');
    rel_name = strrep(rel_name, 'Wam', 'Xa2');
    rel_name = strrep(rel_name, 'ARb', 'Xa5');
    rel_name = strrep(rel_name, 'ARa', 'Xa6');
    rel_name = strrep(rel_name, 'ARf', 'Xa7');
    
    % Peak detectors
    rel_name = strrep(rel_name, 'GC', '');   % REMOVE IF YOU WANT TO LOOK AT QRS DETECTORS.
    rel_name = strrep(rel_name, 'IMS', '');   % REMOVE IF YOU WANT TO LOOK AT PEAK DETECTORS.
    rel_name = strrep(rel_name, 'IMS', 'PDt1');
    rel_name = strrep(rel_name, 'COr', 'PDt2');
    rel_name = strrep(rel_name, 'DCl', 'PDt3');
    
    % Features
    rel_name = strrep(rel_name, 'am', 'Xb2');
    rel_name = strrep(rel_name, 'bwm', 'Xb4');
    rel_name = strrep(rel_name, 'bw', 'Xb1');
    rel_name = strrep(rel_name, 'fm', 'Xb3');
    rel_name = strrep(rel_name, 'pk', 'Xb5');
    rel_name = strrep(rel_name, 'on', 'Xb6');
    rel_name = strrep(rel_name, 'qrsW', 'Xb7');
    rel_name = strrep(rel_name, 'qrsA', 'Xb8');
    rel_name = strrep(rel_name, 'pca', 'Xb9');
    rel_name = strrep(rel_name, 'pulW', 'Xb10');
    
    % get rid of resampling annotations
    rel_name = strrep(rel_name, 'RScubB', '');
    rel_name = strrep(rel_name, 'RSbrgB', '');
    rel_name = strrep(rel_name, 'RSlinB', '');
    rel_name = strrep(rel_name, 'RSlin', '');
    rel_name = strrep(rel_name, 'RScub', '');
    rel_name = strrep(rel_name, 'RSbrg', '');
    % Resampling
    rel_name = strrep(rel_name, 'RSlin', 'RSa1');
    rel_name = strrep(rel_name, 'RScub', 'RSa2');
    rel_name = strrep(rel_name, 'RSbrg', 'RSa3');
    rel_name = strrep(rel_name, 'RScubB', 'RSa4');
    rel_name = strrep(rel_name, 'RSbrgB', 'RSa5');
    
    % RR estimation
    rel_name = strrep(rel_name, 'FTS', 'Ef1');
    rel_name = strrep(rel_name, 'ARS', 'Ef2');
    rel_name = strrep(rel_name, 'ARM', 'Ef3');
    rel_name = strrep(rel_name, 'ARPz', 'Ef5');
    rel_name = strrep(rel_name, 'ARP', 'Ef4');
    rel_name = strrep(rel_name, 'ACF', 'Ef6');
    rel_name = strrep(rel_name, 'WCH', 'Ef7');
    rel_name = strrep(rel_name, 'PKS', 'Et1');
    rel_name = strrep(rel_name, 'ZeX', 'Et2');
    rel_name = strrep(rel_name, 'PZX', 'Et3');
    rel_name = strrep(rel_name, 'CtO', 'Et4');
    rel_name = strrep(rel_name, 'CtA', 'Et5');
    
    % RR Fusion
    rel_name = strrep(rel_name, 'TFu', 'Ft1');
    rel_name = strrep(rel_name, 'PMC', 'Fm2');
    rel_name = strrep(rel_name, 'SFu', 'Fm1');
    rel_name = strrep(rel_name, 'PRC', 'Fm3');
    rel_name = strrep(rel_name, 'SPA', 'Fm4');
    
    % Find signal name
    if ~isempty(strfind(rel_name, 'ekg'))
        alg_names.sigs{alg_no} = 'ECG';
    elseif ~isempty(strfind(rel_name, 'ppg'))
        alg_names.sigs{alg_no} = 'PPG';
    end
    
    % Tidy up
    rel_name = strrep(rel_name, 'ekg_', '');
    rel_name = strrep(rel_name, 'ppg_', '');
    rel_name = strrep(rel_name, '_', ',');
    rel_name = strrep(rel_name, ',,,,,', ',');
    rel_name = strrep(rel_name, ',,,,', ',');
    rel_name = strrep(rel_name, ',,,', ',');
    rel_name = strrep(rel_name, ',,', ',');
    
    if strcmp(rel_name(1),',')
        rel_name = rel_name(2:end);
    end
    
    alg_names.names{alg_no} = rel_name;
    
    clear rel_name
    
end


%% Find elements of algorithms
comp_str = {'Xb', 'Xa', 'RSa', 'PDt', 'Ef', 'Et', 'Fm', 'Ft'};
for comp_no = 1 : length(comp_str)
    eval(['alg_names.meths.' lower(comp_str{comp_no}) ' = nan(length(alg_names.names),1);']);
end
for alg_no = 1 : length(alg_names.names)
    for comp_str = {'Xb', 'Xa', 'RSa', 'PDt', 'Ef', 'Et', 'Fm', 'Ft'}
        rel_el = strfind(alg_names.names{alg_no}, comp_str{1,1});
        if rel_el
            start_el = rel_el+2;
            temp = start_el + 1;
            if temp > length(alg_names.names{alg_no}) || strcmp(alg_names.names{alg_no}(temp), ',')
                temp = temp-1;
            end
            end_el =temp; clear temp
            temp = str2num(alg_names.names{alg_no}(start_el:end_el));
            eval(['alg_names.meths.' lower(comp_str{1,1}), '(alg_no) = temp;']);
        end
        clear rel_el temp start_el end_el
    end
end
clear alg_no comp_str comp_no

%% Save
save(savepath, save_name)

end

function create_table_win_data(up)

%% save setup
save_name = up.paths.filenames.win_data;
savepath = [up.paths.data_save_folder, up.paths.filenames.win_data, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

fprintf('\n--- Creating a Table of Results for each window ');

%% Load algorithm names
loadpath = [up.paths.data_save_folder, up.paths.filenames.alg_names, '.mat'];
load(loadpath); clear loadpath

%% Find out which group each subject is in
subjs = up.paramSet.subj_list;
rel_els = find_rel_subjs(up);
groups = fieldnames(rel_els);
for group_no = 1 : length(groups)
    eval([groups{group_no} '_log = false(length(subjs),1);']);
    eval([groups{group_no} '_log(rel_els.' groups{group_no} ') = true;']);
end

%% Cycle through each algorithm, subject and window
ecg_log = logical(strcmp(alg_names.sigs, 'ECG'));        % is it ecg or ppg
comp_str = {'Xb', 'Xa', 'RSa', 'PDt', 'Ef', 'Et', 'Fm', 'Ft'};

% create storage for raw data
[data.subj, data.ecg_log, data.alg_no] = deal(single([]));
for group_no = 1 : length(groups)
    eval(['data.' groups{group_no} '_log = single([]);']);
end
for comp_no = 1 : length(comp_str)
    eval(['data.m_' lower(comp_str{comp_no}) ' = single([]);']);
end
[data.win_no, data.est, data.ref, data.sqi, data.hr, data.hr_rr, data.snr_log] = deal(single([]));
data.alg_names = cell(0);

% load ref RRs and est RRs for all algorithms for all subjects
counter_no = 0;
for subj = subjs
    
    % Load ref RRs
    loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.ref_rrs, '.mat'];
    load(loadpath);
    
    % Load est RRs
    loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.rrEsts, '.mat'];
    est_RRs = load(loadpath);
    
    % Load SQIs
    loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.sqi, '.mat'];
    load(loadpath);
    
    % Store data
    ppg_hr_rr_ratio = sqi.ppg.hr./rr_ref.v;
    ekg_hr_rr_ratio = sqi.ekg.hr./rr_ref.v;
    for s = 1 : length(alg_names.old_names)
        eval(['rel_data = est_RRs.' alg_names.old_names{s} '.v(:);']);
        start_el = counter_no+1;
        counter_no = counter_no + length(rel_data);
        data.est(start_el:counter_no,1) = rel_data;
        data.alg_no(start_el:counter_no,1) = s;
        if ecg_log(s)
            data.sqi(start_el:counter_no,1) = sqi.ekg.v;
            data.hr(start_el:counter_no,1) = sqi.ekg.hr;
            data.hr_rr(start_el:counter_no,1) = ekg_hr_rr_ratio;
        else
            data.sqi(start_el:counter_no,1) = sqi.ppg.v;
            data.hr(start_el:counter_no,1) = sqi.ppg.hr;
            data.hr_rr(start_el:counter_no,1) = ppg_hr_rr_ratio;
        end
    end
    % variables which are already in the order of the windows
    data.ref = [data.ref; repmat(rr_ref.v(:), [length(alg_names.old_names),1])];
    data.snr_log = [data.snr_log; repmat(rr_ref.snr_log(:), [length(alg_names.old_names),1])];
    win_nos = 1 : length(rr_ref.v);
    data.win_no = [data.win_no; repmat(win_nos(:), [length(alg_names.old_names),1])];
    % variables which are in the order of the algorithms , and need to be changed
    temp = repmat(alg_names.names(:)', [length(rr_ref.v), 1]);
    temp = temp(1:numel(temp)); temp = temp(:);
    data.alg_names = [data.alg_names; temp];
    temp = repmat(ecg_log(:)', [length(rr_ref.v), 1]);
    temp = temp(1:numel(temp)); temp = temp(:);
    data.ecg_log = [data.ecg_log; temp];
    for comp_no = 1 : length(comp_str)
        eval(['temp = repmat(alg_names.meths.' lower(comp_str{comp_no}) '(:)'', [length(rr_ref.v), 1]);']);
        temp = temp(1:numel(temp)); temp = temp(:);
        eval(['data.m_' lower(comp_str{comp_no}) ' = [data.m_' lower(comp_str{comp_no}) '; temp];']);
    end
    % variables which are fixed for this subj
    data.subj = [data.subj; subj*ones(length(rr_ref.v)*length(alg_names.old_names),1)];
    %data.young_log = [data.young_log; young_log(subj)*ones(length(rr_ref.v)*length(alg_names.old_names),1)];
    for group_no = 1 : length(groups)
        eval(['data.' groups{group_no} '_log = [data.' groups{group_no} '_log; ' groups{group_no} '_log(subj)*ones(length(rr_ref.v)*length(alg_names.old_names),1)];']);
    end
end
clear win_nos comp_no s start_el counter_no est_RRs temp rr_ref loadpath subj ecg_log rel_data ppg_hr_rr_ratio sqi ekg_hr_rr_ratio

%data.young_log = logical(data.young_log);
for group_no = 1 : length(groups)
    eval(['data.' groups{group_no} '_log = logical(data.' groups{group_no} '_log);']);
end
data.ecg_log = logical(data.ecg_log);
data.sqi = logical(data.sqi);
data.snr_log = logical(data.snr_log);
data.comb_log = logical(data.sqi & data.snr_log);

%% Save
algorithm_names = data.alg_names;
data = rmfield(data, 'alg_names');
win_data = data;
save(savepath, save_name, 'algorithm_names')

end

function calc_study_and_sub_group_stats(up)

save_name = 'BA_results';
savepath = [up.paths.data_save_folder, up.paths.filenames.global_BA, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

fprintf('\n--- Performing BA analyses for entire cohort and sub-groups ');

%% Load data from entire study
load_name = up.paths.filenames.win_data;
loadpath = [up.paths.data_save_folder, up.paths.filenames.win_data, '.mat'];
load(loadpath, load_name)

%% Identify groups of subjects for which to perform BA analyses:
% i) Entire cohort
% ii) Subgroups
% iii) High and Low Quality
% iv) HR (quintiles)
% v) RR (quintiles)
% vi) HR:RR (quintiles)

group_names = unique(up.paramSet.groups);
groups.entire = true(length(win_data.subj),1);
for group_no = 1 : length(group_names)
    eval(['groups.' group_names{group_no} ' = win_data.' group_names{group_no} '_log;']);
end
%groups.low_qual = ~win_data.sqi;

%% remove these groups for the moment
include_extra_groups = 0;

if include_extra_groups
    no_quantiles = 20;
    win_data.rr = win_data.ref;
    for param = {'hr', 'rr', 'hr_rr'}
        eval(['rel_data = win_data.' param{1,1} ';']);
        for quintile_no = 1 : no_quantiles
            lower_boundary = quantile(rel_data, (quintile_no-1)/no_quantiles);
            if quintile_no == 1
                lower_boundary = lower_boundary - 0.1;
            end
            upper_boundary = quantile(rel_data, quintile_no/no_quantiles);
            eval(['groups.' param{1,1} '_' num2str(quintile_no) ' = win_data.' param{1,1} '>lower_boundary & win_data.' param{1,1} ' <= upper_boundary;']);
            clear upper_boundary lower_boundary
        end
        clear rel_data quintile_no
    end
    clear param
end


%% Find errors
win_data.errors = win_data.est - win_data.ref;

%% BA Analyses

% compute for each group
group_names = fieldnames(groups);
% setup
alg_nos = unique(win_data.alg_no);
no_algs = length(alg_nos);

no_subjs = length(unique(win_data.subj));
for group_no = 1 : length(group_names)
    
    group_name = group_names{group_no,1};
    eval(['rel_group_els = groups.' group_name ';']);
    
    [BA.bias.val, BA.bias.lci, BA.bias.uci, ...
        BA.prec.val, BA.prec.lci, BA.prec.uci, ...
        BA.lloa.lci, BA.lloa.val, BA.lloa.uci, ...
        BA.uloa.val, BA.uloa.lci, BA.uloa.uci, ...
        BA.prop] = deal(nan(no_algs,1));
    
    % for each algorithm ...
    
    for alg_no_no = 1 : no_algs
        
        alg_no = alg_nos(alg_no_no);
        
        % identify relevant data
        if strcmp(group_name, 'low_qual')
            rel_els = rel_group_els & win_data.alg_no == alg_no & ~win_data.sqi & win_data.snr_log;
        else
            % otherwise exclude low quality windows
            rel_els = rel_group_els & win_data.alg_no == alg_no & win_data.sqi & win_data.snr_log;
        end
        
        rel_errors = win_data.errors(rel_els);
        rel_subjs = win_data.subj(rel_els);
        subjs = unique(rel_subjs);
        no_subjs = length(subjs);
        
        % calculate prop of wins
        BA.prop(alg_no) = sum(~isnan(rel_errors)/length(rel_errors));
        
        % calculate bias
        BA.bias.val(alg_no) = nanmean(rel_errors);
        
        % calculate bias CI
        n = no_subjs;    % sample size (no subjects)
        doff = n - 1;               % Degrees of Freedom
        tval = tinv(0.975,doff);    % value from T distribution
        std_diffs = std(rel_errors(~isnan(rel_errors)));    % standard dev of diffs
        std_err_diffs = sqrt((std_diffs^2)/n);      % standard err of diffs
        BA.bias.lci(alg_no) = BA.bias.val(alg_no) - (tval*std_err_diffs);
        BA.bias.uci(alg_no) = BA.bias.val(alg_no) + (tval*std_err_diffs);
        
        % calculate loa
        [p,tbl] = anova1(rel_errors,rel_subjs,'off'); % one-way anova of the differences
        MS_subj = tbl{2,4};
        MS_resid = tbl{3,4};
        no_wins = nan(no_subjs,1);
        for subj_no = 1:no_subjs
            no_wins(subj_no) = sum(rel_subjs == subjs(subj_no));
        end
        denominator = ( (sum(no_wins)^2) - sum(no_wins.^2) ) / ( (n-1)*sum(no_wins) );
        heterogeneity_variance = (MS_subj - MS_resid)/denominator;
        total_variance = heterogeneity_variance + MS_resid;
        BA.lloa.val(alg_no) = BA.bias.val(alg_no) - (1.96*sqrt(total_variance));
        BA.uloa.val(alg_no) = BA.bias.val(alg_no) + (1.96*sqrt(total_variance));
        
        % calculate loa CI
        std_err_lim = sqrt((3*(std_diffs^2))/n);      % standard err of the limit of diffs
        BA.lloa.lci(alg_no) = BA.lloa.val(alg_no) - (tval*std_err_lim);
        BA.lloa.uci(alg_no) = BA.lloa.val(alg_no) + (tval*std_err_lim);
        BA.uloa.lci(alg_no) = BA.uloa.val(alg_no) - (tval*std_err_lim);
        BA.uloa.uci(alg_no) = BA.uloa.val(alg_no) + (tval*std_err_lim);
        
        % calculate precision
        BA.prec.val(alg_no) = BA.uloa.val(alg_no) - BA.lloa.val(alg_no);
        BA.prec.lci(alg_no) = BA.prec.val(alg_no) - (2*std_err_lim*tval);
        BA.prec.uci(alg_no) = BA.prec.val(alg_no) + (2*std_err_lim*tval);
        
        % calculate 2SD
        BA.two_sd.val(alg_no) = 2*sqrt(total_variance);
        
    end
    
    eval(['BA_results.' group_name '= BA;']);
end

%% Save
save(savepath, save_name);

end

function specific_synth_stats(up)

fprintf('\n--- Calculating Statistics for Synthetic Data');

% check to see if it's been done before
save_name = 'synth_results';
savepath = [up.paths.data_save_folder, up.paths.filenames.synth_results, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

% load win_data
load_name = up.paths.filenames.win_data;
loadpath = [up.paths.data_save_folder, up.paths.filenames.win_data, '.mat'];
load(loadpath, load_name)

% find out the rrs, hrs and mods for each subj (from the original data)
load([up.paths.data_load_folder, up.paths.data_load_filename]);
hr = extractfield(data, 'hr');
rr = extractfield(data, 'rr');
mod = extractfield(data, 'group');

% additional data
subjs = unique(win_data.subj);
algs = unique(win_data.alg_no);
errors = win_data.est - win_data.ref;

% find out which modulations can be expected to be extracted by each
% algorithm
load_name = up.paths.filenames.alg_names;
loadpath = [up.paths.data_save_folder, up.paths.filenames.alg_names, '.mat'];
load(loadpath, load_name)
[alg_mods.bw, alg_mods.am, alg_mods.fm] = deal(zeros(length(algs),1));
for alg_no = 1 : length(algs)
    % skip if this is a modulation fusion algorithm
    if ~isnan(alg_names.meths.fm(alg_no)) % || ~isnan(alg_names.meths.tfu(alg_no))
        [alg_mods.bw(alg_no), alg_mods.am(alg_no), alg_mods.fm(alg_no)]= deal(nan);
        continue
    end
    %
    comps.fme = alg_names.meths.xb(alg_no);
    comps.flt = alg_names.meths.xa(alg_no);
    [alg_mods.bw(alg_no), alg_mods.am(alg_no), alg_mods.fm(alg_no)] = identify_alg_mod(comps);
end
% 
% % cycle through each algorithm
% alg_acc = zeros(length(algs),1);
% mod_types = {'bw', 'am', 'fm'};
% for mod_type_no = 1:length(mod_types)
%     mod_type = mod_types{mod_type_no};
%     eval([mod_type '_subjs = find(strcmp(mod, mod_type));']);
% end
% alg_perf = nan(length(algs),1);
% for alg_no = 1 : length(algs)
%     % ignore if this is a fusion algorithm
%     mod_logs = [alg_mods.bw(alg_no), alg_mods.am(alg_no), alg_mods.fm(alg_no)];
%     if sum(isnan(mod_logs)) == length(mod_logs)
%         continue
%     end
%     
%     % cycle through each subject
%     temp = find(mod_logs);
%     no_acc_subjs = zeros(length(temp),1);
%     for temp_no = 1 : length(temp)
%         rel_mod = mod_types(temp(temp_no));
%         eval(['rel_subjs = ' rel_mod{1,1} '_subjs;']);
%         for subj_no = 1 : length(rel_subjs)
%             subj = rel_subjs(subj_no);
%             rel_els = win_data.subj == subj & win_data.alg_no == alg_no;
%             rel_errors = abs(errors(rel_els));
%             no_bad_errors = sum(rel_errors > 1);
%             if no_bad_errors == 0
%                 no_acc_subjs(temp_no) = no_acc_subjs(temp_no) +1;
%             end
%         end
%     end
%     prop_errors = max(no_acc_subjs)/length(rel_subjs);
%     alg_perf(alg_no) = prop_errors;
%     clear no_acc_subjs
% end
% 
% [act_alg_names.names,a,b] = unique(alg_names.names);
% act_alg_names.old_names = alg_names.old_names(a);
% act_alg_names.meths.ef = alg_names.meths.ef(a);
% act_alg_names.meths.et = alg_names.meths.et(a);
% act_alg_names.meths.ft = alg_names.meths.ft(a);
% act_alg_names.meths.xb = alg_names.meths.xb(a);
% a_alg_perf = nan(length(act_alg_names.names),1);
% for aa_no = 1 : length(act_alg_names.names)
%     rel_els = find(strcmp(alg_names.names, act_alg_names.names(aa_no)));
%     temp = alg_perf(rel_els);
%     if sum(isnan(temp)) == length(temp)
%         continue
%     end
%     temp = temp(~isnan(temp));
%     a_alg_perf(aa_no) = max(temp);
% end
% 
% thresh = 0.5;
% synth_results.total_algs = length(a_alg_perf);                              % all algorithms (note that the signal type is irrelevant)
% synth_results.mod_algs = sum(isnan(a_alg_perf));                            % the number of fusion algorithms
% synth_results.remaining_algs = sum(~isnan(a_alg_perf));                     % the number of algorithms after removing fusion algorithms
% synth_results.good_algs = find(~isnan(a_alg_perf) & a_alg_perf >= thresh);  % algorithms deemed to be accurate
% synth_results.bad_algs = find(~isnan(a_alg_perf) & a_alg_perf < thresh);    % algorithms deemed to be inaccurate
% exc_ef_no = 5;
% synth_results.est_9_algs = find(act_alg_names.meths.ef==exc_ef_no);
% synth_results.est_9_bad_algs = find(~isnan(a_alg_perf) & a_alg_perf < thresh & act_alg_names.meths.ef==exc_ef_no);
% synth_results.remaining_aft_est_9_bad_algs = find(~isnan(a_alg_perf) & a_alg_perf < thresh & act_alg_names.meths.ef~=exc_ef_no);
% synth_results.tfu_remaining_aft_est_9_bad_algs = find(~isnan(a_alg_perf) & a_alg_perf < thresh & act_alg_names.meths.ef~=exc_ef_no & ~isnan(act_alg_names.meths.ft) );
% synth_results.remaining_aft_est_9_and_tfu = find(~isnan(a_alg_perf) & a_alg_perf < thresh & act_alg_names.meths.ef~=exc_ef_no & isnan(act_alg_names.meths.ft) );
% synth_results.unique_alg_names = act_alg_names.names;
% synth_results.analysed_algs = find(act_alg_names.meths.ef~=exc_ef_no);
% synth_results.analysed_algs_and_ecg_only = find((act_alg_names.meths.xb==7 | act_alg_names.meths.xb==8) & act_alg_names.meths.ef~=exc_ef_no);
% 
% %% New analysis to get rid of Ef3 and Ef4
% 
% thresh = 0.75;
% mod_synth_results.total_algs = length(a_alg_perf);
% mod_synth_results.mod_algs = sum(isnan(a_alg_perf));
% mod_synth_results.remaining_algs = sum(~isnan(a_alg_perf));
% mod_synth_results.good_algs = find(~isnan(a_alg_perf) & a_alg_perf >= thresh);
% mod_synth_results.bad_algs = find(~isnan(a_alg_perf) & a_alg_perf < thresh);
% mod_synth_results.est_f234_algs = find(~isnan(a_alg_perf) & act_alg_names.meths.ef==4 | act_alg_names.meths.ef==3 | act_alg_names.meths.ef==2);
% mod_synth_results.est_f234_bad_algs = find(~isnan(a_alg_perf) & a_alg_perf < thresh & (act_alg_names.meths.ef==4 | act_alg_names.meths.ef==3 | act_alg_names.meths.ef==2));
% mod_synth_results.remaining_aft_est_f234_bad_algs = find(~isnan(a_alg_perf) & a_alg_perf < thresh & act_alg_names.meths.ef~=4 & act_alg_names.meths.ef~=3 & act_alg_names.meths.ef~=2);
% mod_synth_results.tfu_remaining_aft_est_f234_bad_algs = find(~isnan(a_alg_perf) & a_alg_perf < thresh & (act_alg_names.meths.ef~=4 & act_alg_names.meths.ef~=3 & act_alg_names.meths.ef~=2) & ~isnan(act_alg_names.meths.ft) );
% mod_synth_results.remaining_aft_est_f234_and_tfu = find(~isnan(a_alg_perf) & a_alg_perf < thresh & (act_alg_names.meths.ef~=4 & act_alg_names.meths.ef~=3 & act_alg_names.meths.ef~=2) & isnan(act_alg_names.meths.ft) );
% mod_synth_results.unique_alg_names = act_alg_names.names;
% mod_synth_results.analysed_algs = find(act_alg_names.meths.ef~=4 & act_alg_names.meths.ef~=3 & act_alg_names.meths.ef~=2);
% mod_synth_results.analysed_algs_and_ecg_only = find((act_alg_names.meths.xb==7 | act_alg_names.meths.xb==8) & (act_alg_names.meths.ef~=2 & act_alg_names.meths.ef~=4 & act_alg_names.meths.ef~=3));







%%% Try again using any modulation

% find the 'subjects' which have each modulation
mod_types = {'bw', 'am', 'fm'};
for mod_type_no = 1:length(mod_types)
    mod_type = mod_types{mod_type_no};
    eval([mod_type '_subjs = find(strcmp(mod, mod_type));']);
end

% identify rel_els for each subj and each alg
alg_rel_els = cell(length(algs),1);
for alg_no = 1 : length(algs)
    alg_rel_els{alg_no} = win_data.alg_no == alg_no;
end
subj_rel_els = cell(length(algs),1);
for subj_no = 1 : length(subjs)
    subj_rel_els{subj_no} = win_data.subj == subj_no;
end

sigs = {'ekg', 'ppg'};
sig_rel_els = cell(length(sigs),1);
for sig_no = 1 : length(sigs)
    rel_sig = sigs{sig_no};
    if strcmp(rel_sig, 'ekg')
        sig_rel_els{sig_no} = win_data.ecg_log;
    else
        sig_rel_els{sig_no} = ~win_data.ecg_log;
    end
end

alg_perf = nan(length(algs),6);
% cycle through each signal
for sig_no = 1:length(sigs)
    rel_sig = sigs{sig_no};
    % cycle through each algorithm
    for alg_no = 1 : length(algs)
        % ignore if this is a fusion algorithm
%         if ~isnan(alg_names.meths.fm(alg_no))
%             alg_perf(alg_no, :) = -1;
%             continue
%         end
        fprintf([num2str(alg_no), ', '])
        % cycle through each subject
        [no_rel_subjs, no_acc_subjs] = deal(zeros(3,1));
        for mod_no = 1 : 3
            rel_mod = mod_types(mod_no);
            eval(['rel_subjs = ' rel_mod{1,1} '_subjs;']);
            curr_no_rel_subjs = 0;
            for subj_no = 1 : length(rel_subjs)
                subj = rel_subjs(subj_no);
                rel_els = subj_rel_els{subj} & alg_rel_els{alg_no} & sig_rel_els{sig_no};
                rel_errors = abs(errors(rel_els));
                no_bad_errors = sum(rel_errors > 1);
                if sum(rel_els)==0
                    continue
                end
                curr_no_rel_subjs = curr_no_rel_subjs+1;
                if no_bad_errors == 0
                    no_acc_subjs(mod_no) = no_acc_subjs(mod_no) +1;
                end
                clear rel_errors no_bad_errors rel_els subj
            end
            no_rel_subjs(mod_no) = curr_no_rel_subjs;
            clear subj_no rel_subjs rel_mod curr_no_rel_subjs
        end
        if sum(no_rel_subjs) > 0
            prop_acc = no_acc_subjs./no_rel_subjs;
            if strcmp(rel_sig, 'ekg')
                alg_perf(alg_no, 1:3) = prop_acc;
            else
                alg_perf(alg_no, 4:6) = prop_acc;
            end
        end
        clear no_acc_subjs no_rel_subjs mod_no no_acc_subjs prop_acc
    end
    clear alg_no rel_sig
end
clear sig_no

% alg_perf(alg_perf<0.5) = nan;
% for line_no = 1: length(alg_perf)
%     uselog(line_no) = sum(~isnan(alg_perf(line_no,:)));
% end
% uselog = uselog(:)>0;

[act_alg_names.names,a,b] = unique(alg_names.names);
act_alg_names.old_names = alg_names.old_names(a);
act_alg_names.meths.ef = alg_names.meths.ef(a);
act_alg_names.meths.et = alg_names.meths.et(a);
act_alg_names.meths.ft = alg_names.meths.ft(a);
act_alg_names.meths.xb = alg_names.meths.xb(a);
[all_alg_perf, ppg_alg_perf, ecg_alg_perf] = deal(nan(length(act_alg_names.names),1));
for aa_no = 1 : length(act_alg_names.names)
    rel_els = find(strcmp(alg_names.names, act_alg_names.names(aa_no)));
    temp = alg_perf(rel_els,:);
    
    % all
    if sum(sum(isnan(temp))) == numel(temp)
        continue
    end
    temp1 = temp(~isnan(temp));
    all_alg_perf(aa_no) = max(temp1);
    % ecg
    temp1 = temp(:,1:3);
    if sum(sum(isnan(temp1))) ~= numel(temp1)
        temp1 = temp1(~isnan(temp1));
        ecg_alg_perf(aa_no) = max(temp1);
    end
    % ppg
    temp1 = temp(:,4:6);
    if sum(sum(isnan(temp1))) ~= numel(temp1)
        temp1 = temp1(~isnan(temp1));
        ppg_alg_perf(aa_no) = max(temp1);
    end
    clear temp temp1 rel_els
end

thresh = 0.5;

synth_results.names = act_alg_names.names;
synth_results.total_algs = length(act_alg_names.meths.xb);
synth_results.ecg_only_algs = find(act_alg_names.meths.xb==8 | act_alg_names.meths.xb==7);
synth_results.ppg_only_algs = find(act_alg_names.meths.xb==10);
synth_results.good_algs = find(all_alg_perf >= thresh);    % algorithms deemed to be accurate
synth_results.bad_algs = find(all_alg_perf < thresh);    % algorithms deemed to be inaccurate
synth_results.bad_algs_wo_tfu = find(all_alg_perf < thresh & isnan(act_alg_names.meths.ft));
synth_results.bad_algs_wo_tfu_or_specific_techs = find(all_alg_perf < thresh & act_alg_names.meths.xb~=10 & act_alg_names.meths.xb==9 & act_alg_names.meths.xb==7 & act_alg_names.meths.ef~=5);
synth_results.xb10 = [sum(act_alg_names.meths.xb==10 & all_alg_perf >= thresh), sum(act_alg_names.meths.xb==10)];
synth_results.xb9 = [sum(act_alg_names.meths.xb==9 & all_alg_perf >= thresh), sum(act_alg_names.meths.xb==9)];
synth_results.xb7 = [sum(act_alg_names.meths.xb==7 & all_alg_perf >= thresh), sum(act_alg_names.meths.xb==7)];
synth_results.ef5 = [sum(act_alg_names.meths.ef==5 & all_alg_perf >= thresh), sum(act_alg_names.meths.ef==5)];

% Exclude Xb10 as it performs poorly:
exc_xb_no = 10;
synth_results.names = act_alg_names.names;
synth_results.total_algs = find(act_alg_names.meths.xb~=exc_xb_no);
synth_results.imp_ecg_only_algs = find((act_alg_names.meths.xb==8 | act_alg_names.meths.xb==7) & act_alg_names.meths.xb~=exc_xb_no);
synth_results.good_algs = find(act_alg_names.meths.xb~=exc_xb_no & all_alg_perf >= thresh);    % algorithms deemed to be accurate
synth_results.bad_algs = find(act_alg_names.meths.xb~=exc_xb_no & all_alg_perf < thresh);    % algorithms deemed to be inaccurate
synth_results.bad_algs_wo_tfu = find(act_alg_names.meths.xb~=exc_xb_no & all_alg_perf < thresh & isnan(act_alg_names.meths.ft));

synth_results.total_algs = length(all_alg_perf);                              % all algorithms (note that the signal type is irrelevant)
synth_results.exc_mod_algs = sum(all_alg_perf == -1);                            % the number of fusion algorithms
synth_results.inc_remaining_algs = sum(all_alg_perf ~= -1);                     % the number of algorithms after removing fusion algorithms
synth_results.good_algs = find(all_alg_perf ~= -1 & all_alg_perf >= thresh);  % algorithms deemed to be accurate
synth_results.bad_algs = find(all_alg_perf ~= -1 & all_alg_perf < thresh);    % algorithms deemed to be inaccurate
exc_ef_no = 5;
exc_xb_no = 10;
synth_results.ef_5_algs = find(act_alg_names.meths.ef==exc_ef_no & act_alg_names.meths.xb~=exc_xb_no);
synth_results.ef_5_bad_algs = find(all_alg_perf ~= -1 & all_alg_perf < thresh & act_alg_names.meths.ef==exc_ef_no  & act_alg_names.meths.xb~=exc_xb_no);
synth_results.xb_10_algs = find(act_alg_names.meths.xb==exc_xb_no);
synth_results.xb_10_bad_algs = find(all_alg_perf ~= -1 & all_alg_perf < thresh & act_alg_names.meths.xb==exc_xb_no);
synth_results.xb_9_algs = find(act_alg_names.meths.xb==9 & act_alg_names.meths.ef~=exc_ef_no);
synth_results.xb_9_bad_algs = find(all_alg_perf ~= -1 & all_alg_perf < thresh & act_alg_names.meths.xb==9 & act_alg_names.meths.ef~=exc_ef_no);
synth_results.xb_7_algs = find(act_alg_names.meths.xb==7 & act_alg_names.meths.ef~=exc_ef_no);
synth_results.xb_7_bad_algs = find(all_alg_perf ~= -1 & all_alg_perf < thresh & act_alg_names.meths.xb==7 & act_alg_names.meths.ef~=exc_ef_no);
synth_results.imp_algs = find(act_alg_names.meths.xb~=exc_xb_no);
synth_results.imp_ecg_only_algs = find((act_alg_names.meths.xb==8 | act_alg_names.meths.xb==7));
synth_results.used_algs = find(act_alg_names.meths.ef~=exc_ef_no & act_alg_names.meths.xb~=exc_xb_no);
synth_results.used_ecg_algs = find(act_alg_names.meths.ef~=exc_ef_no & act_alg_names.meths.xb~=exc_xb_no);
synth_results.used_ppg_algs = find(act_alg_names.meths.ef~=exc_ef_no & act_alg_names.meths.xb~=exc_xb_no & act_alg_names.meths.xb~=8 & act_alg_names.meths.xb~=7);

synth_results.remaining_aft_exc_algs = find(all_alg_perf ~= -1 & all_alg_perf < thresh & act_alg_names.meths.ef~=exc_ef_no & act_alg_names.meths.xb~=exc_xb_no);
synth_results.remaining_aft_exc_algs_and_tfu = find(all_alg_perf ~= -1 & all_alg_perf < thresh & act_alg_names.meths.ef~=exc_ef_no & act_alg_names.meths.xb~=exc_xb_no & isnan(act_alg_names.meths.ft) );
synth_results.act_alg_names = act_alg_names;
synth_results.analysed_algs = find(act_alg_names.meths.ef~=exc_ef_no);
synth_results.analysed_algs_and_ecg_only = find((act_alg_names.meths.xb==7 | act_alg_names.meths.xb==8) & act_alg_names.meths.ef~=exc_ef_no);
synth_results.a_alg_perf = all_alg_perf;
synth_results.alg_perf = alg_perf;
synth_results.alg_names = alg_names;

% To look at effects of HR and RR

% for sig_no = 1:length(sigs)
%     % cycle through each algorithm
%     no_errors = nan(length(subjs), length(algs));
%     for alg_no = 1 : length(algs)
%         fprintf([num2str(alg_no), ', '])
%         % cycle through each subject
%         [no_rel_subjs, no_acc_subjs] = deal(zeros(3,1));
%         for subj_no = 1 : length(subjs)
%             subj = subjs(subj_no);
%             rel_els = subj_rel_els{subj} & alg_rel_els{alg_no} & sig_rel_els{sig_no};
%             rel_errors = abs(errors(rel_els));
%             no_errors(subj_no, alg_no) = sum(rel_errors > 1);
%             clear rel_errors rel_els subj
%         end
%         clear subj_no
%     end
%     eval([sigs{sig_no} '_no_errors = no_errors;'])
%     clear alg_no no_errors
% end
% clear sig_no
% 
% rel_cols = alg_names.meths.xb==10;
% rel_rows = 59:87;   % 1-29, 88-122 : BW;   30-58, 123-157 : AM;   59-87, 158-192 : FM;
% plot(mean(ppg_no_errors(rel_rows, rel_cols), 2)), hold on


%% Save
save(savepath, save_name);

end

%% Merged period scripts

function calc_stats_for_rest_and_rec(up)

fprintf('\n--- Calculating Statistics for Combined Rest and Rec Periods ');

%% Copy renamed algorithms
copy_renamed_algorithms(up);

%% Create table of each algorithm, each subject, and each window
create_table_win_data_merged(up);

%% Calculate stats for entire study, and for each sub-group analysis
calc_study_and_sub_group_stats(up);

%% Find specific stats for vortal paper;
specific_vortal_stats(up);
specific_vortal_plots(up);

%% Find traditional stats (e.g. MAE):
calc_trad_stats(up);

end

function copy_renamed_algorithms(up)

%% save setup
save_name = up.paths.filenames.alg_names;
savepath = [up.paths.data_save_folder, up.paths.filenames.alg_names, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

fprintf('\n--- Copying Renamed Algorithms ');

% Specify the file to copy
temp = strrep(up.paths.data_save_folder, '_REST_AND_REC_TEMP', '');
data_load_folder = strrep(temp, '_rest_and_rec', '_rest');
load([data_load_folder, up.paths.filenames.alg_names, '.mat']);

% copy file
save(savepath, save_name)

end

function create_table_win_data_merged_imp(up)

%% save setup
save_name = up.paths.filenames.win_data;
savepath = [up.paths.data_save_folder, up.paths.filenames.win_data_imp, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

fprintf('\n--- Creating a Table of Results for each window ');

%% Load two sets of window data
% rest
root_data_folder = [up.paths.root_folder, 'VORTAL_rest', up.paths.equipment_type, up.paths.slash_direction];
data_save_folder = [root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
load_name = up.paths.filenames.win_data;
loadpath = [data_save_folder, up.paths.filenames.win_data_imp, '.mat'];
rest.data = load(loadpath, load_name);
rest.alg_names = load(loadpath, 'algorithm_names');
rest.data.win_data.rest_log = true(length(rest.data.win_data.subj),1);
% rec
root_data_folder = [up.paths.root_folder, 'VORTAL', up.paths.equipment_type, '_REC', up.paths.slash_direction];
data_save_folder = [root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
load_name = up.paths.filenames.win_data;
loadpath = [data_save_folder, up.paths.filenames.win_data_imp, '.mat'];
rec.data = load(loadpath, load_name); rec.data.win_data.EHV_log = false(length(rec.data.win_data.subj),1);
rec.alg_names = load(loadpath, 'algorithm_names');
rec.data.win_data.subj = translate_subjs_from_rec_to_rest(rec.data.win_data.subj, rest.data);
rec.data.win_data.rest_log = false(length(rec.data.win_data.subj),1);

%% Perform temporal filtering if "REST_AND_REC_TEMP"
max_no_wins = 10;
if strfind(up.paths.data_save_folder, '_TEMP')
    
    % create storage
    params = fieldnames(rest.data.win_data); params = params(:)';
    for param = params
        eval(['new_w_data.' param{1,1} ' = [];']);
    end
    
    % cycle through the specified algorithms
    for alg_no = 1 : length(rel_alg_nos)
        rel_alg_no = rel_alg_nos(alg_no);
        rel_alg_name{alg_no} = rest.alg_names.algorithm_names(min(find(rest.data.win_data.alg_no == rel_alg_no)));
        
        
        % perform for each period
        for period = {'rest', 'rec'}
            eval(['rel_data = ' period{1,1} '.data.win_data;']);
            
            % cycle through each subject
            subjects = unique(rel_data.subj); subjects = subjects(:)';
            for subj = subjects
                
                % identify data corresponding to the relevant algorithm for this subject
                rel_els = rel_data.subj(:) == subj;
                subj_data.est = rel_data.est(rel_els);
                subj_data.ref = rel_data.ref(rel_els);
                subj_data.win_no = rel_data.win_no(rel_els);
                subj_data.sqi = rel_data.sqi(rel_els);
                subj_data.snr_log = rel_data.snr_log(rel_els);
                
                % cycle through the number of wins to include in the filtering
                different_no_wins = 1:max_no_wins;
                for no_wins = different_no_wins
                    % filter data
                    [new_est, new_ref] = med_filt_data(subj_data, no_wins);
                    temp = ones(length(new_est),1);
                    % set all parameters to be the first value
                    counter_no = length(new_w_data.subj);
                    for param = params
                        eval(['temp2 = rel_data.' param{1,1} '(rel_els);']);
                        temp2 = temp2(1);
                        eval(['new_w_data.' param{1,1} '(counter_no+1:counter_no+length(temp),1) = temp2*temp;']);
                    end
                    % reset those parameters which vary within a subject
                    new_w_data.alg_no(counter_no+1:counter_no+length(temp),1) = ones(length(temp),1);
                    new_w_data.est(counter_no+1:counter_no+length(temp),1) = new_est;
                    new_w_data.ref(counter_no+1:counter_no+length(temp),1) = new_ref;
                    new_w_data.snr_log(counter_no+1:counter_no+length(temp),1) = subj_data.snr_log;
                    new_w_data.win_no(counter_no+1:counter_no+length(temp),1) = subj_data.win_no;
                    new_w_data.sqi(counter_no+1:counter_no+length(temp),1) = subj_data.sqi;
                    
                end
                
            end
            
        end
        
    end
    
    new_w_data.sqi = logical(new_w_data.sqi);
    new_w_data.comb_log = logical(new_w_data.sqi & new_w_data.snr_log);
    
    win_data = new_w_data;
    algorithm_names = rel_alg_name;
    
    %% Re-make alg_names file appropriately:
    save_name = up.paths.filenames.alg_names;
    savepath = [up.paths.data_save_folder, up.paths.filenames.alg_names, '.mat'];
    
    clear alg_names
    alg_names.meths = struct;
    counter_no = 1;
    for no_wins = different_no_wins
        alg_names.names{counter_no,1} = [rel_alg_name{alg_no}{1,1}, '_' num2str(no_wins)];
        alg_names.sigs{counter_no,1} = rel_alg_sigs{alg_no};
        alg_names.no_wins(counter_no,1) = no_wins;
        counter_no = counter_no+1;
    end
    alg_names.old_names = alg_names.names;
    
    save(savepath, save_name)
    
else
    
    %% Combine into one win_data file
    params = fieldnames(rest.data.win_data);
    for param = params(:)'
        eval(['win_data.' param{1,1} ' = [rest.data.win_data.' param{1,1} '(:); rec.data.win_data.' param{1,1} '(:)];'])
    end
    win_data.comb_log = logical(win_data.sqi & win_data.snr_log);
    
    % merge algorithm names
    algorithm_names = rest.alg_names.algorithm_names;
    algorithm_names(end+1: end+length(rec.alg_names.algorithm_names)) = rec.alg_names.algorithm_names;
    
end

%% Save
save_name = up.paths.filenames.win_data;
savepath = [up.paths.data_save_folder, up.paths.filenames.win_data_imp, '.mat'];
save(savepath, save_name, 'algorithm_names')

end

function create_table_win_data_merged(up)

%% save setup
save_name = up.paths.filenames.win_data;
savepath = [up.paths.data_save_folder, up.paths.filenames.win_data, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

fprintf('\n--- Creating a Table of Results for each window ');

%% Load algorithm names
loadpath = [up.paths.data_save_folder, up.paths.filenames.alg_names, '.mat'];
load(loadpath); clear loadpath

%% Load two sets of window data
% rest
root_data_folder = [up.paths.root_folder, 'VORTAL', up.paths.equipment_type, '_rest', up.paths.slash_direction];
data_save_folder = [root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
load_name = up.paths.filenames.win_data;
loadpath = [data_save_folder, up.paths.filenames.win_data, '.mat'];
rest.data = load(loadpath, load_name);
rest.alg_names = load(loadpath, 'algorithm_names');
rest.data.win_data.rest_log = true(length(rest.data.win_data.subj),1);
% rec
root_data_folder = [up.paths.root_folder, 'VORTAL', up.paths.equipment_type, '_REC', up.paths.slash_direction];
data_save_folder = [root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
load_name = up.paths.filenames.win_data;
loadpath = [data_save_folder, up.paths.filenames.win_data, '.mat'];
rec.data = load(loadpath, load_name); rec.data.win_data.EHV_log = false(length(rec.data.win_data.subj),1);
rec.alg_names = load(loadpath, 'algorithm_names');
rec.data.win_data.subj = translate_subjs_from_rec_to_rest(rec.data.win_data.subj, rest.data);
rec.data.win_data.rest_log = false(length(rec.data.win_data.subj),1);

%% Perform temporal filtering if "REST_AND_REC_TEMP"
max_no_wins = 16;
if strfind(up.paths.data_save_folder, '_TEMP')
    
    % create storage
    params = fieldnames(rest.data.win_data); params = params(:)';
    for param = params
        eval(['new_w_data.' param{1,1} ' = [];']);
    end
    
    % specify which algorithms to perform this on:
    rel_alg_nos = [203,212];   % ECG Est7,MFu2, PPG Est7,MFu2  (small bias)
    %rel_alg_nos = [206, 214];  % ECG Est6,MFu2, PPG Est4,MFu2  (large bias)
    rel_alg_sigs = {'ECG'; 'PPG'};
    
    % cycle through the specified algorithms
    for alg_no = 1 : length(rel_alg_nos)
        rel_alg_no = rel_alg_nos(alg_no);
        rel_alg_name{alg_no} = rest.alg_names.algorithm_names(min(find(rest.data.win_data.alg_no == rel_alg_no)));
        
        % perform for each period
        for period = {'rest', 'rec'}
            eval(['rel_data = ' period{1,1} '.data.win_data;']);
            
            % cycle through each subject
            subjects = unique(rel_data.subj); subjects = subjects(:)';
            for subj = subjects
                
                % identify data corresponding to the relevant algorithm for this subject
                rel_els = rel_data.alg_no == rel_alg_no & rel_data.subj(:) == subj;
                subj_data.est = rel_data.est(rel_els);
                subj_data.ref = rel_data.ref(rel_els);
                subj_data.win_no = rel_data.win_no(rel_els);
                subj_data.sqi = rel_data.sqi(rel_els);
                subj_data.snr_log = rel_data.snr_log(rel_els);
                
                % cycle through the number of wins to include in the filtering
                different_no_wins = 1:max_no_wins;
                for no_wins = different_no_wins
                    % filter data
                    [new_est, new_ref] = med_filt_data(subj_data, no_wins);
                    temp = ones(length(new_est),1);
                    % set all parameters to be the first value
                    counter_no = length(new_w_data.subj);
                    for param = params
                        eval(['temp2 = rel_data.' param{1,1} '(rel_els);']);
                        temp2 = temp2(1);
                        eval(['new_w_data.' param{1,1} '(counter_no+1:counter_no+length(temp),1) = temp2*temp;']);
                    end
                    % reset those parameters which vary within a subject
                    new_w_data.alg_no(counter_no+1:counter_no+length(temp),1) = (max_no_wins*(alg_no-1))+no_wins*ones(length(temp),1);
                    new_w_data.est(counter_no+1:counter_no+length(temp),1) = new_est;
                    new_w_data.ref(counter_no+1:counter_no+length(temp),1) = new_ref;
                    new_w_data.win_no(counter_no+1:counter_no+length(temp),1) = subj_data.win_no;
                    new_w_data.sqi(counter_no+1:counter_no+length(temp),1) = subj_data.sqi;
                    new_w_data.snr_log(counter_no+1:counter_no+length(temp),1) = subj_data.snr_log;
                    
                end
                
            end
            
        end
        
    end
    
    new_w_data.sqi = logical(new_w_data.sqi);
    new_w_data.snr_log = logical(new_w_data.snr_log);
    
    win_data = new_w_data;
    win_data.comb_log = logical(win_data.sqi & win_data.snr_log);
    algorithm_names = rel_alg_name;
    
    %% Re-make alg_names file appropriately:
    save_name = up.paths.filenames.alg_names;
    savepath = [up.paths.data_save_folder, up.paths.filenames.alg_names, '.mat'];
    
    clear alg_names
    alg_names.meths = struct;
    counter_no = 1;
    for alg_no = 1 : length(rel_alg_nos)
        for no_wins = different_no_wins
            alg_names.names{counter_no,1} = [rel_alg_name{alg_no}{1,1}, '_' num2str(no_wins)];
            alg_names.sigs{counter_no,1} = rel_alg_sigs{alg_no};
            alg_names.no_wins(counter_no,1) = no_wins;
            counter_no = counter_no+1;
        end
    end
    alg_names.old_names = alg_names.names;
    
    save(savepath, save_name)
    
else
    
    %% Combine into one win_data file
    params = fieldnames(rest.data.win_data);
    for param = params(:)'
        eval(['win_data.' param{1,1} ' = [rest.data.win_data.' param{1,1} '(:); rec.data.win_data.' param{1,1} '(:)];'])
    end
    % remove the alg_no field because it hasn't accounted for the numbering
    % being different between rest and rec:
    %win_data = rmfield(win_data, 'alg_no');
    
    win_data.comb_log = logical(win_data.sqi & win_data.snr_log);
    
    % merge algorithm names
    algorithm_names = rest.alg_names.algorithm_names;
    algorithm_names(end+1: end+length(rec.alg_names.algorithm_names)) = rec.alg_names.algorithm_names;
    
    %temp = strcat(algorithm_names, num2str(win_data.ecg_log));
    % insert a new alg_no field
    %[alg_names,b,new_alg_nos] = unique(temp);
    %win_data.alg_no = single(new_alg_nos);
    
end

%% Save
save_name = up.paths.filenames.win_data;
savepath = [up.paths.data_save_folder, up.paths.filenames.win_data, '.mat'];
save(savepath, save_name, 'algorithm_names')

end

function [new_est, new_ref] = med_filt_data(subj_data, no_wins)

% cycle through each window, and find med filtered est and ref values
win_nos = unique(subj_data.win_no); win_nos = win_nos(:)';
[new_ref, new_est] = deal(nan(length(win_nos),1));
for win_no = win_nos
    start_win_no = win_no;
    end_win_no = win_no + no_wins - 1;
    rel_els = subj_data.win_no >= start_win_no & subj_data.win_no <= end_win_no ...
        & subj_data.sqi & ~isnan(subj_data.ref) & subj_data.snr_log;
    new_est(win_no) = nanmedian(subj_data.est(rel_els));
    new_ref(win_no) = nanmedian(subj_data.ref(rel_els));
end
end

function calc_trad_stats(up)

save_name = 'trad_stats';
savepath = [up.paths.data_save_folder, up.paths.filenames.study_stats, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

fprintf('\n--- Calculating Traditional Stats');

%% Load alg names
load_name = up.paths.filenames.alg_names;
loadpath = [up.paths.data_save_folder, up.paths.filenames.alg_names, '.mat'];
load(loadpath, load_name);

%% Load windata
load_name = up.paths.filenames.win_data;
loadpath = [up.paths.data_save_folder, up.paths.filenames.win_data, '.mat'];
load(loadpath, load_name)

%% Find stats for each algorithm, for each group
subjs = up.paramSet.subj_list;
rel_els = find_rel_subjs(up);
groups = fieldnames(rel_els);
groups{end+1} = 'entire';
%groups{end+1} = 'low_qual';                                                % uncomment to look at the low quality windows

for group_no = 1 : length(groups)
    
    % identify rel rows
    if ~strcmp(groups{group_no}, 'entire') & ~strcmp(groups{group_no}, 'low_qual')
        eval(['rel_group_rows = win_data.' groups{group_no} '_log;']);
    else
        rel_group_rows = true(length(win_data.subj),1);
    end
    
    for alg_no = 1 : length(alg_names.names)
        
        alg_rows = rel_group_rows & win_data.alg_no == alg_no;
        if ~strcmp(groups{group_no}, 'low_qual')
            rel_rows = alg_rows & win_data.sqi & win_data.snr_log;
        else
            rel_rows = alg_rows & ~win_data.sqi & win_data.snr_log;
        end
        
        % find errors
        error = win_data.est(rel_rows) - win_data.ref(rel_rows);
        
        % Mean Differences
        mean_diff = nanmean(error);
        
        % SD of Differences
        rel_els = ~isnan(error);
        std_diff = std(error(rel_els));
        
        % MAE
        abs_error = abs(error);
        mae = nanmean(abs_error);
        rel_els = ~isnan(abs_error);
        sdae = std(abs_error(rel_els));
        
        % RMSE
        rmse = sqrt(nanmean(error.^2));
        
        % Percentage Error
        mean_ref_rr = nanmean(win_data.ref(rel_rows));
        percerr = 100*nanmean(abs_error/mean_ref_rr);
        
        % Prop Wins ref rr
        
        total_wins_alg = sum(alg_rows);
        total_wins_good_sqi = sum(alg_rows & win_data.sqi);
        total_wins_est = sum(alg_rows & ~isnan(win_data.est));
        total_wins_ref = sum(alg_rows & ~isnan(win_data.ref) & win_data.snr_log);
        total_wins_good_sqi_and_ref = sum(alg_rows & ~isnan(win_data.ref) & win_data.snr_log & win_data.sqi);
        total_wins_all = sum(alg_rows & ~isnan(win_data.est) & ~isnan(win_data.ref) & win_data.comb_log);
        
        % prop of wins which had a reference which had a good sqi
        prop_wins_good_sqi = sum(alg_rows & ~isnan(win_data.ref) & win_data.comb_log)/sum(alg_rows & ~isnan(win_data.ref) & win_data.snr_log);
        % prop of wins which had a reference and a good sqi which had an estimate
        prop_wins_est = sum(alg_rows & ~isnan(win_data.est) & ~isnan(win_data.ref) & win_data.comb_log)/sum(alg_rows & ~isnan(win_data.ref) & win_data.comb_log);
        % prop of wins which had a reference
        prop_wins_ref = total_wins_ref/total_wins_alg;
        % prop of wins which had a reference and a good sqi
        prop_wins_good_sqi_and_ref = total_wins_good_sqi_and_ref/total_wins_alg;
        % prop of wins included in the analysis
        prop_wins_inc = sum(alg_rows & ~isnan(win_data.ref) & win_data.comb_log)/total_wins_alg;
        prop_wins_all = total_wins_all/total_wins_alg;
        
        % Prop with abs error of < 1 bpm
        prop_acc = sum(abs_error < 1)/sum(~isnan(abs_error));
        
        % store each statistic
        for stat = {'percerr', 'mae', 'sdae', 'rmse', 'prop_wins_good_sqi', 'prop_wins_ref', 'prop_wins_est', 'prop_wins_all', 'prop_wins_good_sqi_and_ref', 'total_wins_alg', 'total_wins_good_sqi', 'total_wins_ref', 'total_wins_est', 'total_wins_est', 'total_wins_all', 'total_wins_good_sqi_and_ref', 'prop_acc', 'mean_diff', 'std_diff', 'mean_ref_rr'}
            eval(['trad_stats.' groups{group_no} '.' stat{1,1} '(alg_no) = ' stat{1,1} ';']);
        end
        
    end
    
end

%% Save
save(savepath, save_name)

end


function [bw, am, fm] = identify_alg_mod(comps)


if ~isnan(comps.fme)
    a = comps.fme;
    % Feature-based
    switch a
        case 1
            bw = 1; am = 0; fm = 0;
        case 2
            bw = 0; am = 1; fm = 0;
        case 3
            bw = 0; am = 0; fm = 1;
        case 4
            bw = 1; am = 0; fm = 0;
        case 5
            bw = 1; am = 1; fm = 0;
        case 6
            bw = 1; am = 1; fm = 0;
        case 7
            bw = 0; am = 0; fm = 1;
        case 8
            bw = 0; am = 1; fm = 1;
        case 9
            bw = 1; am = 1; fm = 1;
        case 10
            bw = 0; am = 0; fm = 1;
    end
    
else
    % Filter-based
    a = comps.flt;
    switch a
        case 1
            bw = 1; am = 0; fm = 0;
        case 2
            bw = 0; am = 1; fm = 0;
        case 3
            bw = 0; am = 0; fm = 1;
        case 4
            bw = 1; am = 1; fm = 1;
    end
end

end

function new_subj_list = translate_subjs_from_rec_to_rest(old_subj_list, rest)

% find out what numbers the YHVs were in rest:
rest_subjs = unique(rest.win_data.subj(rest.win_data.young_log));

% load rec subjects
curr_rec_subjs = unique(old_subj_list);

for subj_no = 1 : length(curr_rec_subjs)
    new_rec_subjs(subj_no) = rest_subjs(curr_rec_subjs(subj_no));
    new_subj_list(old_subj_list == curr_rec_subjs(subj_no)) = new_rec_subjs(subj_no);
end

end

function create_table_win_data_merged_raw_clin(up)

%% save setup
save_name = up.paths.filenames.win_data_raw_clin;
savepath = [up.paths.data_save_folder, up.paths.filenames.win_data_raw_clin, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

fprintf('\n--- Merging Tables of Results for Raw and Clinical data');

%% Load algorithm names
% Clin
% loadpath = [up.paths.data_save_folder, up.paths.filenames.alg_names, '.mat'];
% clin = load(loadpath); clear loadpath
% Raw
%raw_save_folder = strrep(up.paths.data_save_folder, '_clin', '');
% loadpath = [raw_save_folder, up.paths.filenames.alg_names, '.mat'];
% raw = load(loadpath); clear loadpath

%% Load four sets of window data
% rest (clin)
root_data_folder = [up.paths.root_folder, 'VORTAL_clin', up.paths.slash_direction];
data_save_folder = [root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
load_name = up.paths.filenames.win_data;
loadpath = [data_save_folder, up.paths.filenames.win_data, '.mat'];
clin.rest.data = load(loadpath, load_name);
clin.rest.alg_names = load(loadpath, 'algorithm_names');
clin.rest.data.win_data.rest_log = true(length(clin.rest.data.win_data.subj),1);
clin.rest.data.win_data.raw_log = false(length(clin.rest.data.win_data.subj),1);
% rec (clin)
root_data_folder = [up.paths.root_folder, 'VORTAL_clin_rec', up.paths.slash_direction];
data_save_folder = [root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
load_name = up.paths.filenames.win_data;
loadpath = [data_save_folder, up.paths.filenames.win_data, '.mat'];
clin.rec.data = load(loadpath, load_name); clin.rec.data.win_data.EHV_log = false(length(clin.rec.data.win_data.subj),1);
clin.rec.alg_names = load(loadpath, 'algorithm_names');
clin.rec.data.win_data.subj = translate_subjs_from_rec_to_rest(clin.rec.data.win_data.subj, clin.rest.data);
clin.rec.data.win_data.rest_log = false(length(clin.rec.data.win_data.subj),1);
clin.rec.data.win_data.raw_log = false(length(clin.rec.data.win_data.subj),1);
% rest (raw)
root_data_folder = [up.paths.root_folder, 'VORTAL', up.paths.slash_direction];
data_save_folder = [root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
load_name = up.paths.filenames.win_data;
loadpath = [data_save_folder, up.paths.filenames.win_data, '.mat'];
raw.rest.data = load(loadpath, load_name);
raw.rest.alg_names = load(loadpath, 'algorithm_names');
raw.rest.data.win_data.rest_log = true(length(raw.rest.data.win_data.subj),1);
raw.rest.data.win_data.raw_log = true(length(raw.rest.data.win_data.subj),1);
% rec (raw)
root_data_folder = [up.paths.root_folder, 'VORTAL_rec', up.paths.slash_direction];
data_save_folder = [root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
load_name = up.paths.filenames.win_data;
loadpath = [data_save_folder, up.paths.filenames.win_data, '.mat'];
raw.rec.data = load(loadpath, load_name); raw.rec.data.win_data.EHV_log = false(length(raw.rec.data.win_data.subj),1);
raw.rec.alg_names = load(loadpath, 'algorithm_names');
raw.rec.data.win_data.subj = translate_subjs_from_rec_to_rest(raw.rec.data.win_data.subj, raw.rest.data);
raw.rec.data.win_data.rest_log = false(length(raw.rec.data.win_data.subj),1);
raw.rec.data.win_data.raw_log = true(length(raw.rec.data.win_data.subj),1);

%% Combine into one win_data file
params = fieldnames(raw.rest.data.win_data);
for param = params(:)'
    eval(['win_data.' param{1,1} ' = [raw.rest.data.win_data.' param{1,1} '(:);' ...
        'raw.rec.data.win_data.' param{1,1} '(:);' ...
        'clin.rest.data.win_data.' param{1,1} '(:);' ...
        'clin.rec.data.win_data.' param{1,1} '(:);];'])
end

win_data.comb_log = logical(win_data.sqi & win_data.snr_log);

% merge algorithm names
algorithm_names = raw.rest.alg_names.algorithm_names;
algorithm_names(end+1: end+length(raw.rec.alg_names.algorithm_names)) = raw.rec.alg_names.algorithm_names;
algorithm_names(end+1: end+length(clin.rest.alg_names.algorithm_names)) = clin.rest.alg_names.algorithm_names;
algorithm_names(end+1: end+length(clin.rec.alg_names.algorithm_names)) = clin.rec.alg_names.algorithm_names;

win_data_raw_clin = win_data;

%% Save
save_name = up.paths.filenames.win_data_raw_clin;
savepath = [up.paths.data_save_folder, up.paths.filenames.win_data_raw_clin, '.mat'];
save(savepath, save_name, 'algorithm_names')

end

function create_table_win_data_merged_raw_clin_imp(up)

%% save setup
save_name = up.paths.filenames.win_data_raw_clin_imp;
savepath = [up.paths.data_save_folder, up.paths.filenames.win_data_raw_clin_imp, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

fprintf('\n--- Merging Tables of Results for Raw and Clinical data');

%% Load four sets of window data
% rest (clin)
root_data_folder = [up.paths.root_folder, 'VORTAL_clin', up.paths.slash_direction];
data_save_folder = [root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
load_name = up.paths.filenames.win_data;
loadpath = [data_save_folder, up.paths.filenames.win_data_imp, '.mat'];
clin.rest.data = load(loadpath, load_name);
clin.rest.alg_names = load(loadpath, 'algorithm_names');
clin.rest.data.win_data.rest_log = true(length(clin.rest.data.win_data.subj),1);
clin.rest.data.win_data.raw_log = false(length(clin.rest.data.win_data.subj),1);
% rec (clin)
root_data_folder = [up.paths.root_folder, 'VORTAL_clin_rec', up.paths.slash_direction];
data_save_folder = [root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
load_name = up.paths.filenames.win_data;
loadpath = [data_save_folder, up.paths.filenames.win_data_imp, '.mat'];
clin.rec.data = load(loadpath, load_name); clin.rec.data.win_data.EHV_log = false(length(clin.rec.data.win_data.subj),1);
clin.rec.alg_names = load(loadpath, 'algorithm_names');
clin.rec.data.win_data.subj = translate_subjs_from_rec_to_rest(clin.rec.data.win_data.subj, clin.rest.data);
clin.rec.data.win_data.rest_log = false(length(clin.rec.data.win_data.subj),1);
clin.rec.data.win_data.raw_log = false(length(clin.rec.data.win_data.subj),1);
% rest (raw)
root_data_folder = [up.paths.root_folder, 'VORTAL', up.paths.slash_direction];
data_save_folder = [root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
load_name = up.paths.filenames.win_data;
loadpath = [data_save_folder, up.paths.filenames.win_data_imp, '.mat'];
raw.rest.data = load(loadpath, load_name);
raw.rest.alg_names = load(loadpath, 'algorithm_names');
raw.rest.data.win_data.rest_log = true(length(raw.rest.data.win_data.subj),1);
raw.rest.data.win_data.raw_log = true(length(raw.rest.data.win_data.subj),1);
% rec (raw)
root_data_folder = [up.paths.root_folder, 'VORTAL_rec', up.paths.slash_direction];
data_save_folder = [root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
load_name = up.paths.filenames.win_data;
loadpath = [data_save_folder, up.paths.filenames.win_data_imp, '.mat'];
raw.rec.data = load(loadpath, load_name); raw.rec.data.win_data.EHV_log = false(length(raw.rec.data.win_data.subj),1);
raw.rec.alg_names = load(loadpath, 'algorithm_names');
raw.rec.data.win_data.subj = translate_subjs_from_rec_to_rest(raw.rec.data.win_data.subj, raw.rest.data);
raw.rec.data.win_data.rest_log = false(length(raw.rec.data.win_data.subj),1);
raw.rec.data.win_data.raw_log = true(length(raw.rec.data.win_data.subj),1);

%% Combine into one win_data file
params = fieldnames(raw.rest.data.win_data);
for param = params(:)'
    eval(['win_data.' param{1,1} ' = [raw.rest.data.win_data.' param{1,1} '(:);' ...
        'raw.rec.data.win_data.' param{1,1} '(:);' ...
        'clin.rest.data.win_data.' param{1,1} '(:);' ...
        'clin.rec.data.win_data.' param{1,1} '(:);];'])
end

win_data.comb_log = logical(win_data.sqi & win_data.snr_log);

% merge algorithm names
algorithm_names = raw.rest.alg_names.algorithm_names;
algorithm_names(end+1: end+length(raw.rec.alg_names.algorithm_names)) = raw.rec.alg_names.algorithm_names;
algorithm_names(end+1: end+length(clin.rest.alg_names.algorithm_names)) = clin.rest.alg_names.algorithm_names;
algorithm_names(end+1: end+length(clin.rec.alg_names.algorithm_names)) = clin.rec.alg_names.algorithm_names;

win_data_raw_clin_imp = win_data;

%% Save
save_name = up.paths.filenames.win_data_raw_clin_imp;
savepath = [up.paths.data_save_folder, up.paths.filenames.win_data_raw_clin_imp, '.mat'];
save(savepath, save_name, 'algorithm_names')

end



%% Unknown scripts

function print_win_table

%% Make table
headers = fieldnames(tab); headers = headers';
table_data = nan(tot_wins, length(headers));
for s = 1 : length(headers)
    eval(['table_data(:,s) = tab.' headers{s} ';']);
end

% Save to csv and mat files
save_name = [up.paths.filenames.win_results_table];
savepath = [up.paths.tables_save_folder, save_name];
csvwrite([savepath, '.csv'],table_data)
save(savepath, 'table_data')

save_name = [up.paths.filenames.win_results_headers];
savepath = [up.paths.tables_save_folder, save_name];
xlswrite(savepath,headers)
save(savepath, 'headers')

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
    for comp_str = {'Sb', 'Sa', 'RSa', 'PDt', 'Ea', 'Eb', 'Fa', 'Fb'}
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

function calc_rep_meas_BA(up)

% Find relevant subjs
rel_els = find_rel_subjs(up);
groups = fieldnames(rel_els);


% Skip if this processing has been done previously
if ~up.analysis.redo_stats
    return
end

for group_no = 1 : length(groups)
    
    group = groups{group_no};
    eval(['rel_group_els = rel_els.' group ';']);
    
    % load global data
    load_name = ['data_', group];
    loadpath = [up.paths.data_save_folder, up.paths.filenames.group_data, group, '.mat'];
    load(loadpath, load_name);
    eval(['comp_rel_data = ' load_name ';']);
    
    % for each algorithm ...
    alg_names = fieldnames(comp_rel_data);
    [BA.bias.val, BA.bias.lci, BA.bias.uci, ...
        BA.prec.val, BA.prec.lci, BA.prec.uci, ...
        BA.lloa.lci, BA.lloa.val, BA.lloa.uci, ...
        BA.uloa.val, BA.uloa.lci, BA.uloa.uci] = deal(nan(length(alg_names),1));
    
    for alg_no = 1 : length(alg_names)
        
        % extract rel data
        eval(['alg_rel_data = comp_rel_data.' alg_names{alg_no} ';']);
        
        % calculate errors (in bpm)
        alg_rel_data.errors = alg_rel_data.est - alg_rel_data.ref;
        
        % calculate bias
        BA.bias.val(alg_no) = nanmean(alg_rel_data.errors);
        
        % calculate bias CI
        n = length(unique(alg_rel_data.subj));    % sample size (no subjects)
        doff = n - 1;               % Degrees of Freedom
        tval = tinv(0.975,doff);    % value from T distribution
        std_diffs = std(alg_rel_data.errors(~isnan(alg_rel_data.errors)));    % standard dev of diffs
        std_err_diffs = sqrt((std_diffs^2)/n);      % standard err of diffs
        BA.bias.lci(alg_no) = BA.bias.val(alg_no) - (tval*std_err_diffs);
        BA.bias.uci(alg_no) = BA.bias.val(alg_no) + (tval*std_err_diffs);
        
        % calculate loa
        [p,tbl] = anova1(alg_rel_data.errors,alg_rel_data.subj,'off'); % one-way anova of the differences
        MS_subj = tbl{2,4};
        MS_resid = tbl{3,4};
        no_wins = nan(length(rel_group_els),1);
        for subj_no = 1:length(rel_group_els)
            no_wins(subj_no) = sum(alg_rel_data.subj == rel_group_els(subj_no));
        end
        denominator = ( (sum(no_wins)^2) - sum(no_wins.^2) ) / ( (n-1)*sum(no_wins) );
        heterogeneity_variance = (MS_subj - MS_resid)/denominator;
        total_variance = heterogeneity_variance + MS_resid;
        BA.lloa.val(alg_no) = BA.bias.val(alg_no) - (1.96*sqrt(total_variance));
        BA.uloa.val(alg_no) = BA.bias.val(alg_no) + (1.96*sqrt(total_variance));
        
        % calculate loa CI
        std_err_lim = sqrt((3*(std_diffs^2))/n);      % standard err of the limit of diffs
        BA.lloa.lci(alg_no) = BA.lloa.val(alg_no) - (tval*std_err_lim);
        BA.lloa.uci(alg_no) = BA.lloa.val(alg_no) + (tval*std_err_lim);
        BA.uloa.lci(alg_no) = BA.uloa.val(alg_no) - (tval*std_err_lim);
        BA.uloa.uci(alg_no) = BA.uloa.val(alg_no) + (tval*std_err_lim);
        
        % calculate precision
        BA.prec.val(alg_no) = BA.uloa.val(alg_no) - BA.lloa.val(alg_no);
        BA.prec.lci(alg_no) = BA.prec.val(alg_no) - (2*std_err_lim*tval);
        BA.prec.uci(alg_no) = BA.prec.val(alg_no) + (2*std_err_lim*tval);
        
        % calculate 2SD
        BA.two_sd.val(alg_no) = 2*sqrt(total_variance);
        
    end
    
    save_name = 'BA';
    savepath = [up.paths.data_save_folder, 'BA_', group, '.mat'];
    % Save Global BA to file
    save_or_append_data
    
end

end

function calc_entire_sample_stats(win_stats_analyses, up)

if ~up.analysis.redo_stats
    return
end

global_stats_analyses = {'accuracy'}; %, 'coverage'};
% cycle through stats analyses
% Find relevant subjs
rel_els = find_rel_subjs(up);
groups = fieldnames(rel_els);

for group_no = 1 : length(groups)    % for each group (YHVs, EHVs, global)
    for global_stats_analysis_no = 1 : length(win_stats_analyses)
        
        % Skip if this processing has been done previously
        save_name = win_stats_analyses{global_stats_analysis_no}(1:3);
        savepath = [up.paths.data_save_folder, up.paths.filenames.group_stats, groups{group_no}, '.mat'];
        if ~up.analysis.redo_stats
            exist_log = check_exists(savepath, save_name);
            if exist_log
                continue
            end
        end
        
        % Load data from entire sample
        global_stats_data = cell(1);
        eval(['temp_rel_els = rel_els.' groups{group_no} ';']);
        rel_subjs = up.paramSet.subj_list(temp_rel_els);
        for subj_no = 1:length(rel_subjs)
            loadpath = [up.paths.data_save_folder, num2str(rel_subjs(subj_no)), up.paths.filenames.stats, '.mat'];
            global_stats_data{subj_no,1} = load(loadpath, ['win_' global_stats_analyses{global_stats_analysis_no}]);
        end
        
        % collate list of RR algorithms
        if ~exist('rr_alg_names', 'var')
            eval(['temp = global_stats_data{1,1}.win_' win_stats_analyses{global_stats_analysis_no} ';'])
            temp2 = fieldnames(temp); temp2 = temp2{1};
            eval(['rr_alg_names = fieldnames(temp.' temp2 ');']);
        end
        
        
        
        %% Med and IQR of Accuracy Stats:
        eval(['acc_stats = fieldnames(global_stats_data{subj_no,1}.win_' win_stats_analyses{global_stats_analysis_no} ');']);
        for s = 1 : length(acc_stats)
            
            acc_stat = acc_stats{s};
            
            % ignore the raw data
            if sum(strcmp(acc_stat, {'error', 'est_rr', 'ref_rr'}))
                continue
            end
            
            % make blank vector
            eval([ '[' win_stats_analyses{global_stats_analysis_no}(1:3) '.' acc_stat '.med, acc.' acc_stat '.iqr] = deal(nan(length(rr_alg_names),1));']);
            
            % Perform each stats analysis on entire sample:
            for rr_alg_no = 1 : length(rr_alg_names)
                
                % extract stat for each subject
                stat_data = nan(length(rel_subjs),1);
                for subj_no = 1:length(rel_subjs)
                    eval(['temp = global_stats_data{subj_no,1}.win_' win_stats_analyses{global_stats_analysis_no} '.' acc_stat '.' rr_alg_names{rr_alg_no} ';'])
                    stat_data(subj_no) = temp;
                end
                % Ignore NaNs   (CHECK: worth thinking about)
                stat_data = stat_data(~isnan(stat_data));
                
                % Calcs + Store results for this algorithm
                eval([ win_stats_analyses{global_stats_analysis_no}(1:3) '.' acc_stat '.med(rr_alg_no) = median(stat_data);']);
                eval([ win_stats_analyses{global_stats_analysis_no}(1:3) '.' acc_stat '.iqr(rr_alg_no) = iqr(stat_data);']);
                eval([ win_stats_analyses{global_stats_analysis_no}(1:3) '.' acc_stat '.lq(rr_alg_no) = quantile(stat_data, 0.25);']);
                eval([ win_stats_analyses{global_stats_analysis_no}(1:3) '.' acc_stat '.uq(rr_alg_no) = quantile(stat_data, 0.75);']);
                
            end
            
        end
        eval([ win_stats_analyses{global_stats_analysis_no}(1:3) '.rr_alg_names = rr_alg_names;']);
        
    end
    %% Save stats to file
    save_or_append_data
end

end

function calc_win_stats(win_stats_analyses, up)

for subj = up.paramSet.subj_list
    
    % cycle through stats analyses
    for win_stats_analysis_no = 1 : length(win_stats_analyses)
        
        %% Skip if this processing has been done previously
        save_name = ['win_' win_stats_analyses{win_stats_analysis_no}];
        savepath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.stats, '.mat'];
        if ~up.analysis.redo_stats
            exist_log = check_exists(savepath, save_name);
            if exist_log
                continue
            end
        end
        
        %% Load estimated RRs for this subj:
        if ~exist('est_RRs', 'var')
            loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.rrEsts, '.mat'];
            est_RRs = load(loadpath);
        end
        
        %% Load reference RRs for this subj:
        if ~exist('rr_ref', 'var')
            loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.ref_rrs, '.mat'];
            load(loadpath); % loads rr_ref
        end
        
        %% Perform each stats analysis on each subj:
        temp = feval(win_stats_analyses{win_stats_analysis_no}, est_RRs, rr_ref, up);
        eval([save_name ' = temp;']);
        
        %% Save stats to file
        save_or_append_data
        
    end
    
    clear rr_ref est_RRs
    
end

end

function stats = accuracy(est_RRs, ref_RRs, up)

estRR_names = fieldnames(est_RRs);

for estRR_no = 1 : length(estRR_names)
    
    eval(['rel_estRRs = est_RRs.' estRR_names{estRR_no} ';']);
    
    % calc stat for each window over which RR was estimated
    error = nan(length(ref_RRs.t),1);
    for win_no = 1 : length(ref_RRs.t)
        error(win_no) = rel_estRRs.v(win_no) - ref_RRs.v(win_no);
    end
    
    % Mean Differences
    mean_diff = nanmean(error);
    
    % SD of Differences
    rel_els = ~isnan(error);
    std_diff = std(error(rel_els));
    
    % Mean ref RR
    mean_ref_rr = nanmean(ref_RRs.v);
    
    % Est RR
    est_rr = rel_estRRs.v;
    % Ref RR
    ref_rr = ref_RRs.v;
    
    % MAE
    abs_error = abs(error);
    mae = nanmean(abs_error);
    rel_els = ~isnan(abs_error);
    sdae = std(abs_error(rel_els));
    
    % RMSE
    rmse = sqrt(nanmean(error.^2));
    
    % Percentage Error
    percerr = 100*nanmean(abs(est_rr-ref_rr)/mean_ref_rr);
    
    % Prop Wins
    prop_wins = sum(~isnan(error))/length(error);
    
    % Prop with abs error of < 1 bpm
    prop_acc = sum(abs_error < 1)/sum(~isnan(abs_error));
    
    % store calc data
    eval(['stats.error.' estRR_names{estRR_no} ' = error;']);
    
    % store each statistic
    for stat = {'percerr', 'mae', 'sdae', 'rmse', 'prop_wins', 'prop_acc', 'mean_diff', 'std_diff', 'mean_ref_rr', 'est_rr', 'ref_rr'}
        eval(['stats.' stat{1,1} '.' estRR_names{estRR_no} ' = ' stat{1,1} ';']);
    end
    
    clear rel_estRRs error
    
end

end
