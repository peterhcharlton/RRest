function fuse_rr(up)
%FUSE_RR fuses RRs estimated using a variety of techniques, as specified
% in PC's literature review.
%	            fuse_rr(up)
%
%	Inputs:
%		data            data files should be stored in the specified format
%       up              universal parameters structure
%
%	Outputs:
%       	for each subject, n:
%       n_rrEsts.m      - a file of RR estimates, with fused estimates
%                           added to those previously estimated.
%

fprintf('\n--- Fusing RRs ');

%% Identify fusable RR ests
possible_rr_ests = find_possible_rr_ests(up);

%% Mod fusion
% Identify relevant rrEsts to fuse for mods:
SFu_rr_ests = find_rr_ests_for_mods(possible_rr_ests);
SPA_rr_ests = find_rr_ests_for_mods(possible_rr_ests, {'WCH', 'FTS', 'ACF', 'ARS', 'ARM'});
PMC_rr_ests = find_rr_ests_for_mods(possible_rr_ests, {'ARP'});
PRC_rr_ests = find_rr_ests_for_mods(possible_rr_ests, {'ARP'});
rel_sub_comps = up.al.sub_components.fus_mod;
for subj = up.paramSet.subj_list
    loaded_this_subj_rr_ests = 0;
    for sub_comp_no = 1 : length(rel_sub_comps)
        % Skip if this processing has been done previously
        eval(['temp = ' rel_sub_comps{sub_comp_no} '_rr_ests;']);
        for fusable_rr_est_no = 1:length(temp.bef)
            part1 = temp.bef{fusable_rr_est_no};
            part2 = temp.aft{fusable_rr_est_no};
            eval(['save_name = ''' part1, '_', part2, '_', rel_sub_comps{sub_comp_no} ''';']);
            savepath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.rrEsts, '.mat'];
            check_exists
            % load data if it hasn't yet been loaded
            if ~loaded_this_subj_rr_ests
                rrEsts = load([up.paths.data_save_folder, num2str(subj), up.paths.filenames.rrEsts]);
                loaded_this_subj_rr_ests = 1;
            end
            % Identify relevant data
            feat_mods = {'bw', 'am', 'fm'};
            filt_mods = {'ARb', 'ARa', 'ARf'};
            rr_est_names = fieldnames(rrEsts);
            for mod_no = 1:3
                mod = feat_mods(mod_no);
                rel_rr_est = [temp.bef{fusable_rr_est_no}, mod{1,1}, temp.aft{fusable_rr_est_no}];
                if sum(strcmp(rr_est_names, rel_rr_est))
                    % then it is a feature mod
                    eval(['rel_data.' mod{1,1} ' = rrEsts.' temp.bef{fusable_rr_est_no}, mod{1,1}, temp.aft{fusable_rr_est_no} ';']);
                    eval(['rel_data.' mod{1,1} ' = rrEsts.' temp.bef{fusable_rr_est_no}, mod{1,1}, temp.aft{fusable_rr_est_no} ';']);
                    eval(['rel_data.' mod{1,1} ' = rrEsts.' temp.bef{fusable_rr_est_no}, mod{1,1}, temp.aft{fusable_rr_est_no} ';']);
                else
                    mod = filt_mods(mod_no);
                    save_mod = feat_mods(mod_no);
                    eval(['rel_data.' save_mod{1,1} ' = rrEsts.' temp.bef{fusable_rr_est_no}, mod{1,1}, temp.aft{fusable_rr_est_no} ';']);
                    eval(['rel_data.' save_mod{1,1} ' = rrEsts.' temp.bef{fusable_rr_est_no}, mod{1,1}, temp.aft{fusable_rr_est_no} ';']);
                    eval(['rel_data.' save_mod{1,1} ' = rrEsts.' temp.bef{fusable_rr_est_no}, mod{1,1}, temp.aft{fusable_rr_est_no} ';']);
                end
            end
            % Do fusion
            temp_rr = feval(rel_sub_comps{sub_comp_no}, rel_data, up); clear rel_data
            % store this series of rrs:
            eval([save_name, ' = temp_rr;']);
            %% Save RRs to file
            save_or_append_data
        end
    end
end

%% Temp fusion
% Identify fusable RR ests a second time since some more will now have been
% made.
if sum(strcmp(up.al.options.fuse_rr, 'fus_temp'))
    possible_rr_ests = find_possible_rr_ests(up);
    % Identify relevant rrEsts to fuse for mods:
    TFu_rr_ests = possible_rr_ests;
    rel_sub_comps = up.al.sub_components.fus_temp;
    for subj = up.paramSet.subj_list
        loaded_this_subj_rr_ests = 0;
        for sub_comp_no = 1 : length(rel_sub_comps)
            % Skip if this processing has been done previously
            eval(['temp = ' rel_sub_comps{sub_comp_no} '_rr_ests;']);
            for fusable_rr_est_no = 1:length(temp)
                eval(['save_name = ''' temp{fusable_rr_est_no}, '_', rel_sub_comps{sub_comp_no} ''';']);
                savepath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.rrEsts, '.mat'];
                check_exists
                % load data if it hasn't yet been loaded
                if ~loaded_this_subj_rr_ests
                    rrEsts = load([up.paths.data_save_folder, num2str(subj), up.paths.filenames.rrEsts]);
                    loaded_this_subj_rr_ests = 1;
                end
                % Identify relevant data
                eval(['rel_data = rrEsts.' temp{fusable_rr_est_no} ';']);
                % Do fusion
                temp_rr = feval(rel_sub_comps{sub_comp_no}, rel_data, up);
                % store this series of rrs:
                eval([save_name, ' = temp_rr;']);
                %% Save RRs to file
                save_or_append_data
            end
        end
        clear rrEsts
    end
end

end

function possible_rr_ests = find_possible_rr_ests(up)

% Extract list of RR ests from the first pt:
subj = up.paramSet.subj_list(1);
loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.rrEsts, '.mat'];
if exist(loadpath, 'file')
    filecontents = whos('-file', loadpath);
    rrEsts = extractfield(filecontents, 'name');
else
    warning(['No respiratory signals found for Subject ', num2str(subj) '.'])
end

% eliminate those ests which certainly cannot be mod fused because they are a temp fusion:
for s = 1 : length(up.al.sub_components.fus_temp)
    rrEsts_to_exclude = strfind(rrEsts, up.al.sub_components.fus_temp{s});
    log_keep = ones(length(rrEsts_to_exclude),1);
    for s = 1 : length(rrEsts_to_exclude)
        if ~isempty(rrEsts_to_exclude{s})
            log_keep(s) = 0;
        end
    end
end
if ~exist('log_keep', 'var')
    log_keep = true(length(rrEsts),1);
end
possible_rr_ests = cell(length(rrEsts),1);
count_no = 0;
for s= 1 : length(rrEsts)
    if log_keep(s)
        count_no = count_no+1;
        possible_rr_ests(count_no) = rrEsts(s);
    end
end
possible_rr_ests = possible_rr_ests(1:count_no);

end
function ind_mod_rr_ests = find_rr_ests_for_mods(new_rrEsts, exc_criteria)

if nargin == 2
    % then exclude the rrEsts according to the criteria in argument 2
    for s = 1 : length(exc_criteria)
        rrEsts_to_exclude = strfind(new_rrEsts, exc_criteria{s});
        log_keep = ones(length(rrEsts_to_exclude),1);
        for s = 1 : length(rrEsts_to_exclude)
            if isempty(rrEsts_to_exclude{s})
                log_keep(s) = 0;
            end
        end
    end
    rr_ests_to_keep = cell(length(new_rrEsts),1);
    count_no = 0;
    for s= 1 : length(new_rrEsts)
        if log_keep(s)
            count_no = count_no+1;
            rr_ests_to_keep(count_no) = new_rrEsts(s);
        end
    end
    rr_ests_to_keep = rr_ests_to_keep(1:count_no);
    new_rrEsts = rr_ests_to_keep;
end

% AM, BW, FM
temp = strfind(new_rrEsts, 'bw_'); % ignore things like "bwm" (instead of "bw") which don't have all three mods.
[ind_mod_rr_ests.bef, ind_mod_rr_ests.aft] = deal(cell(0));
for s = 1 : length(temp)
    if ~isempty(temp{s})
        curr_ind = length(ind_mod_rr_ests.bef)+1;
        orig_temp = new_rrEsts(s);
        if iscell(orig_temp)
            orig_temp = orig_temp{1};
        end
        ind_mod_rr_ests.bef{curr_ind} = orig_temp(1:(temp{s}-1));
        ind_mod_rr_ests.aft{curr_ind} = orig_temp((temp{s}+2):end);
    end
end

% ARb, ARa, ARf
temp = strfind(new_rrEsts, 'ARb_'); % ignore things like "bwm" (instead of "bw") which don't have all three mods.
%[ind_mod_rr_ests.bef, ind_mod_rr_ests.aft] = deal(cell(0));
for s = 1 : length(temp)
    if ~isempty(temp{s})
        curr_ind = length(ind_mod_rr_ests.bef)+1;
        orig_temp = new_rrEsts(s);
        if iscell(orig_temp)
            orig_temp = orig_temp{1};
        end
        ind_mod_rr_ests.bef{curr_ind} = orig_temp(1:(temp{s}-1));
        ind_mod_rr_ests.aft{curr_ind} = orig_temp((temp{s}+3):end);
    end
end

end