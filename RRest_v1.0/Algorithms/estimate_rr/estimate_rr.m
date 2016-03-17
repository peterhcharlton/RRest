function estimate_rr(up)
%ESTIMATE_RR estimates RR from respiratory signals using each possible set
% of RR estimation options, as specified in PC's literature review.
%	            estimate_rr(up)
%
%	Inputs:
%		data            data files should be stored in the specified format
%       up              universal parameters structure
%
%	Outputs:
%       	for each subject, n:
%       n_rrEsts.m      - a file of RR estimates
%

fprintf('\n--- Estimating RRs ');
%% Extract list of resp signals from the first pt:
subj = up.paramSet.subj_list(1);
loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.respSigs];
if exist(loadpath, 'file')
    filecontents = whos('-file', loadpath);
    respSigs = extractfield(filecontents, 'name');
else
    warning(['No respiratory signals found for Subject ', num2str(subj) '.'])
end

%% Cycle through each patient
for subj = up.paramSet.subj_list
    loaded_this_subj_respSigs = 0;
    %% Make window timings if necessary
    identify_subj_wins(subj, up);
    %% Cycle through each resp signal
    for respSig_no = 1:length(respSigs)
        for option_no = 1 : length(up.al.options.estimate_rr)
            % Skip if this processing has been done previously
            save_name = [ respSigs{respSig_no} '_' up.al.options.estimate_rr{option_no} ];
            savepath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.rrEsts, '.mat'];
            exist_log = check_exists(savepath, save_name);
            if exist_log
                continue
            end
            
            % load data if it hasn't yet been loaded
            if ~loaded_this_subj_respSigs
                % Signals
                load([up.paths.data_save_folder, num2str(subj), up.paths.filenames.respSigs]);
                % Window timings
                load([up.paths.data_save_folder, num2str(subj), up.paths.filenames.win_timings, '.mat']);
                loaded_this_subj_respSigs = 1;
            end
            % Identify the relevant respSig data
            eval(['rel_data = ' respSigs{respSig_no} ';']);
            %% Calculate RR from this resp sig using each option for estimating RR
            if (length(rel_data.t) == 1 && isnan(rel_data.t)) || sum(isnan(rel_data.v))==length(rel_data.v)
                temp_rr.t = mean([wins.t_start(:)' ; wins.t_end(:)']); temp_rr.t = temp_rr.t(:);
                temp_rr.v = nan(length(temp_rr.t),1);
            else
                temp_rr = feval(up.al.options.estimate_rr{option_no}, rel_data, wins, up);
            end
            % store this series of rrs:
            eval([respSigs{respSig_no} '_' up.al.options.estimate_rr{option_no} ' = temp_rr;']);
            clear temp_rr
            %% Save RRs to file
            save_or_append_data
        end
    end
    clear ekg* ppg* wins            % clear resp sigs
end

end