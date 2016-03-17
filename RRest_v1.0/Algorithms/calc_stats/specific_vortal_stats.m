function specific_vortal_stats(up)

save_name = 'vortal_results';
savepath = [up.paths.data_save_folder, up.paths.filenames.vortal_results, '.mat'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        continue
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

%% Find out which group each subject is in
subjs = up.paramSet.subj_list;
rel_subj_els = find_rel_subjs(up);
groups = fieldnames(rel_subj_els);
young_log = false(length(subjs),1);
for subj_no = 1 : length(subjs)
    if sum(rel_subj_els.young == subj_no)
        young_log(subj_no) = true;
    end
end

%% YHV SQI results
load_name = up.paths.filenames.win_data;
loadpath = [up.paths.data_save_folder, up.paths.filenames.win_data, '.mat'];
load(loadpath, load_name);
ecg_alg_no = win_data.alg_no(min(find(win_data.ecg_log == 1)));
ppg_alg_no = win_data.alg_no(min(find(win_data.ecg_log == 0)));
rel_paw_els = win_data.alg_no == ecg_alg_no & win_data.ecg_log == 1;
rel_ppg_els = win_data.alg_no == ppg_alg_no & win_data.ecg_log == 0;
rel_ekg_els = win_data.alg_no == ecg_alg_no & win_data.ecg_log == 1;

for sig = {'paw', 'ekg', 'ppg'}
    for subj_no = 1 : length(rel_subj_els.young)
        eval(['rel_els = rel_' sig{1,1} '_els;']);
        subj = rel_subj_els.young(subj_no);
        rel_sqi_els = find(rel_els & win_data.subj == subj);
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

%% YHV BA results
rel_res = BA_results.young;

% Top ten algorithms for each of ECG and PPG
no_algs_to_disp = 10;
if no_algs_to_disp > length(find(strcmp(alg_names.sigs, 'ECG')))
    no_algs_to_disp = length(find(strcmp(alg_names.sigs, 'ECG')));
end
for sig = {'ekg', 'ppg'}
    if strcmp(sig{1,1}, 'ekg')
        rel_alg_nos = find(strcmp(alg_names.sigs, 'ECG'));
    else
        rel_alg_nos = find(strcmp(alg_names.sigs, 'PPG'));
    end
    [~, order] = sort(rel_res.prec.val(rel_alg_nos));
    
    eval(['top_ten.' sig{1,1} '. alg_nos = rel_alg_nos(order(1:no_algs_to_disp));']);
    for param = {'prec', 'bias', 'lloa', 'uloa'}
        for param2 = {'val', 'lci', 'uci'}
            eval(['top_ten.' sig{1,1} '.' param{1,1} '.' param2{1,1} ' = rel_res.' param{1,1} '.' param2{1,1} '(top_ten.' sig{1,1} '.alg_nos);']);
        end
    end
    eval(['top_ten.' sig{1,1} '.prop = 100*rel_res.prop(top_ten.' sig{1,1} '.alg_nos);']);
end

vortal_results.top_ten = top_ten;


%% Find top algorithm for each technique
% want to find the most precise algorithm for each of the techniques: find
% it's new name and precision.
alg_component = fieldnames(alg_names.meths);
rel_BA_res = BA_results.young;

for sig = {'ECG', 'PPG'}
    [tech_res.prec, tech_res.lci, tech_res.uci] = deal([]);
    [tech_res.meth_name, tech_res.alg_name] = deal(cell(0));
    counter_no = 0;
    for alg_comp_no = 1 : length(alg_component)
        
        rel_alg_component = alg_component{alg_comp_no};
        eval(['rel_meth_vals = alg_names.meths.' rel_alg_component ';']);
        meths = unique(rel_meth_vals); meths = meths(~isnan(meths));
        for meth_no = 1 : length(meths)
            meth = meths(meth_no);
            rel_els = find(rel_meth_vals == meth & strcmp(alg_names.sigs, sig{1,1}));
            precs = rel_BA_res.prec.val(rel_els);
            [prec, el] = min(precs);
            lci = rel_BA_res.prec.lci(rel_els(el));
            uci = rel_BA_res.prec.uci(rel_els(el));
            alg_el = rel_els(el);
            counter_no = counter_no+1;
            if ~isempty(alg_el)
                tech_res.alg_name{counter_no,1} = alg_names.names{alg_el};
                tech_res.meth_name{counter_no,1} = [rel_alg_component, num2str(meth)];
                tech_res.prec(counter_no,1) = prec;
                tech_res.lci(counter_no,1) = lci;
                tech_res.uci(counter_no,1) = uci;
                tech_res.res{counter_no,1} = [num2str(prec,4), ' (' num2str(lci,4), '-', num2str(uci,4) ')'];
            end
        end
    end
    eval([sig{1,1}, '_tech_res = tech_res;']);
    clear tech_res
end

vortal_results.tech_res.ekg = ECG_tech_res;
vortal_results.tech_res.ppg = PPG_tech_res;


%% Save
save(savepath, save_name);

end