function specific_vortal_stats(up)

save_name = 'vortal_results';
savepath = [up.paths.data_save_folder, up.paths.filenames.vortal_results, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

fprintf('\n--- Performing algorithm ranking analyses');

%% Load BA data
load_name = 'BA_results';
loadpath = [up.paths.data_save_folder, up.paths.filenames.global_BA, '.mat'];
load(loadpath, load_name);

%% Load alg names
load_name = up.paths.filenames.alg_names;
loadpath = [up.paths.data_save_folder, up.paths.filenames.alg_names, '.mat'];
load(loadpath, load_name);

subjs = up.paramSet.subj_list;
%% Find out which group each subject is in (now redundant)
rel_subj_els = find_rel_subjs(up, 'group');
groups = fieldnames(rel_subj_els);
young_log = false(length(subjs),1);
for subj_no = 1 : length(subjs)
    if sum(rel_subj_els.young == subj_no)
        young_log(subj_no) = true;
    end
end

%% SQI results
load_name = up.paths.filenames.win_data;
loadpath = [up.paths.data_save_folder, up.paths.filenames.win_data, '.mat'];
load(loadpath, load_name);
ecg_alg_no = win_data.alg_no(min(find(win_data.ecg_log == 1)));
ppg_alg_no = win_data.alg_no(min(find(win_data.ecg_log == 0)));

sigs = [up.paramSet.ekg_sigs(:)', up.paramSet.ppg_sigs(:)', 'paw'];
for sig = sigs
    if ~strcmp(sig{1,1}(1:3), 'paw')
        eval(['rel_sig_log = win_data.' sig{1,1} '_log;']);
        eval(['rel_alg_no = min(win_data.alg_no(win_data.' sig{1,1} '_log));']);
        rel_log = rel_sig_log & win_data.alg_no == rel_alg_no;
    else
        eval(['rel_sig_log = win_data.' up.paramSet.ekg_sigs{1,1} '_log;']);
        eval(['rel_alg_no = min(win_data.alg_no(win_data.' up.paramSet.ekg_sigs{1,1} '_log));']);
        rel_log = rel_sig_log & win_data.alg_no == rel_alg_no;
    end
    
    for subj_no = 1 : length(subjs)
        subj = subjs(subj_no);
        rel_sqi_els = find(rel_log & win_data.subj == subj);
        if ~strcmp(sig, 'paw')
            rel_sqi = win_data.sqi(rel_sqi_els) & win_data.snr_log(rel_sqi_els);
        else
            rel_sqi = win_data.snr_log(rel_sqi_els);
        end
        eval([sig{1,1} '_sqi.high.raw(subj_no) = sum(rel_sqi);']);
        eval([sig{1,1} '_sqi.low.raw(subj_no) = length(rel_sqi) - sum(rel_sqi);']);
        eval([sig{1,1} '_sqi.all.raw(subj_no) = length(rel_sqi);']);
    end
    for cat = {'high', 'low', 'all'}
        eval(['temp_data = ' sig{1,1} '_sqi.' cat{1,1} '.raw;']);
        eval([sig{1,1} '_sqi.' cat{1,1} '.med = median(temp_data);']);
        eval([sig{1,1} '_sqi.' cat{1,1} '.lq = quantile(temp_data, 0.25);']);
        eval([sig{1,1} '_sqi.' cat{1,1} '.uq = quantile(temp_data, 0.75);']);
    end
    eval(['vortal_results.' sig{1,1} '_sqi = ' sig{1,1} '_sqi;']);
end

%% Save
save(savepath, save_name);

end