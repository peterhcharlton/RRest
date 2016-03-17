function RS(option, up)
%RS resamples feat-based respiratory signals at a regular sampling rate.
%	            RS(option, up)
%
%	Inputs:
%		option          the option which has led to this function being used
%       up              universal parameters structure
%

fprintf('\n--- Resampling ');
log_int_respSig = 1;             % Has value 1 unless this is a final respiratory signal

for subj = up.paramSet.subj_list
    
    %% Skip if this processing has been done previously
    iden_resp_sig_file_ending
    savepath = [up.paths.data_save_folder, num2str(subj), ending];
    filecontents = whos('-file', savepath);
    var_names = extractfield(filecontents, 'name');
    rel_log = zeros(size(var_names));
    for s = 1 : length(var_names)
        if strfind(var_names{s}, [option(1:3), up.paths.filenames.feat_meas])
            rel_log(s) = 1;
        end
    end
    rel_var_names = var_names(logical(rel_log));
    for rel_var_name_no = 1 : length(rel_var_names)
        for current_opt_no = 1 : length(up.al.options.RS)
            eval(['save_name = ''' option(1:3), up.paths.filenames.resampler, up.al.options.RS{current_opt_no}, rel_var_names{rel_var_name_no}(4:end) ''';']);
            exist_log = check_exists(savepath, save_name);
            if exist_log
                continue
            end
            
            %% Load relevant data
            loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.int_respSigs];
            rel_name = rel_var_names{rel_var_name_no};
            load(loadpath, rel_name);
            eval(['fs_orig = ' rel_name, '.fs;']);
            eval(['feat_data = ' rel_name, ';']);
            
            %% Resample each waveform using each option
            fs = up.paramSet.resample_fs;
            if ~strcmp(up.al.options.RS{current_opt_no}(end), 'B')
                % resample waveform
                resampled_data = feval(up.al.options.RS{current_opt_no}, feat_data, fs, fs_orig);
                resampled_data.fs = fs;
            else
                resampled_data = feval(up.al.options.RS{current_opt_no}(1:(end-1)), feat_data, fs, fs_orig);
                resampled_data.fs = fs;
                % BPF if specified by option
                try
                    resampled_data = bpf_signal_to_remove_non_resp_freqs(resampled_data, resampled_data.fs, up);
                catch
                    % if there aren't enough samples, then don't BPF:
                    resampled_data = resampled_data;
                end
            end
            % store data
            eval([save_name ' = resampled_data;']); clear resampled_data
            
            %% Save processed data
            save_or_append_data
        end
    end
    
end

end

function resampled_data = cub(feat_data, fs, fs_orig)

if length(feat_data.t) <2
    resampled_data.t = nan;
    resampled_data.v = nan;
    return
end

resampled_data.t = feat_data.t(1):(1/fs):feat_data.t(end);
resampled_data.v = interp1(feat_data.t, feat_data.v, resampled_data.t, 'pchip');

end

function resampled_data = lin(feat_data, fs, fs_orig)

if length(feat_data.t) <2
    resampled_data.t = nan;
    resampled_data.v = nan;
    return
end

resampled_data.t = feat_data.t(1):(1/fs):feat_data.t(end);
resampled_data.v = interp1(feat_data.t, feat_data.v, resampled_data.t, 'linear');

end

function resampled_data = brg(feat_data, fs, fs_orig)
% Using the algorithm described in:
% R. D. Berger et al., “An efficient algorithm for spectral analysis of heart rate variability.,” IEEE Trans. Biomed. Eng., vol. 33, no. 9, pp. 900–4, Sep. 1986.

% Sampling rate for resampled output is fs (this is currently 5, but
% perhaps it should be more like 2 Hz, according to Berger).

% Identify timings of local windows
dummy = feat_data.t(1):(1/fs):feat_data.t(end);
win_centres = dummy(2:(end-1));
win_starts = dummy(1:(end-2));
win_ends = dummy(3:end);

%% Find instantaneous value at each window timing
resampled_data.t = win_centres;
resampled_data.v = nan(length(win_centres),1);

inst_vals = nan(length(win_centres),1);
for win_no = 1 : length(win_centres)
    % Find timing of QRS immediately before the window:
    temp = feat_data.t - win_starts(win_no); temp = temp(temp<=0);
    dummy = -1*min(abs(temp));
    if isempty(dummy)
        resampled_data.v(win_no) = nan;
        continue
    end
    bef_qrs.t = win_starts(win_no) + dummy;
    bef_qrs.i = find(feat_data.t == bef_qrs.t);
    % Find timing of QRS immediately after the window:
    temp = feat_data.t - win_ends(win_no); temp = temp(temp>=0);
    dummy = min(abs(temp));
    if isempty(dummy)
        resampled_data.v(win_no) = nan;
        continue
    end
    aft_qrs.t = win_ends(win_no) + dummy;
    aft_qrs.i = find(feat_data.t == aft_qrs.t);
    % Find QRS's contained within this window:
    inside_qrs.i = find(feat_data.t > win_starts(win_no) & feat_data.t < win_ends(win_no));
    inside_qrs.t = feat_data.t(inside_qrs.i);
    % for each heart beat find instantaneous value
    sections.start_t = [win_starts(win_no); inside_qrs.t];
    sections.end_t = [inside_qrs.t; win_ends(win_no)];
    [section_vals, no_rr_ints] = deal(nan(length(sections.start_t),1));
    for section_no = 1 : length(sections.start_t)
        % section duration
        section_length = sections.end_t(section_no) - sections.start_t(section_no);
        % find qrs immediately before this section
        temp = feat_data.t - sections.start_t(section_no); temp = temp(temp<=0);
        dummy = -1*min(abs(temp));
        rel_qrs_bef.t = sections.start_t(section_no) + dummy;
        rel_qrs_bef.i = find(0.001*round(1000*feat_data.t) == 0.001*round(1000*rel_qrs_bef.t));   % to correct for decimal pt errors
        % Find value at this qrs time
        rel_val = feat_data.v(rel_qrs_bef.i);
        % find qrs immediately after this section
        temp = feat_data.t - sections.start_t(section_no); temp = temp(temp>0);
        dummy = min(abs(temp));
        rel_qrs_aft.t = sections.start_t(section_no) + dummy;
        rel_qrs_aft.i = find(feat_data.t == rel_qrs_aft.t);
        % RR interval during this section
        rr_int = rel_qrs_aft.t - rel_qrs_bef.t;
        % No of RR intervals in this section
        no_rr_ints(section_no) = section_length/rr_int;
        % Find val for this section
        section_vals(section_no) = rel_val;     % This is slightly different to the algorithm to allow for BW and AM
    end
    resampled_data.v(win_no) = (fs/2)*sum(section_vals.*no_rr_ints);
end

end