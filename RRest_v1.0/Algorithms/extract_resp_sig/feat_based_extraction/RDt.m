function RDt(option, up)
%RDt detects QRS waves in ECG signals.
%	            RDt(option, up)
%
%	Inputs:
%		option          the option which has led to this function being used
%       up              universal parameters structure
%
%	Outputs:
%       ...
%

fprintf('\n--- Detecting QRS waves ');
log_int_respSig = 1;             % Has value 1 unless this is a final respiratory signal

for subj = up.paramSet.subj_list
    %% Skip if this processing has been done previously
    iden_resp_sig_file_ending
    savepath = [up.paths.data_save_folder, num2str(subj), ending];
    filecontents = whos('-file', savepath);
    var_names = extractfield(filecontents, 'name');
    temp = strfind(var_names, [option(1:3), up.paths.filenames.elim_vhf]); temp = temp{:};
    rel_var_names = var_names(temp);
    for rel_var_name_no = 1 : length(rel_var_names)
        for current_opt_no = 1 : length(up.al.options.RDt)
            eval(['save_name = ''' option(1:3), up.paths.filenames.qrss, up.al.options.RDt{current_opt_no}, rel_var_names{rel_var_name_no}(4:end) ''';']);
            exist_log = check_exists(savepath, save_name);
            if exist_log
                continue
            end
            
            %% Load relevant data
            loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.int_respSigs];
            rel_name = rel_var_names{rel_var_name_no};
            load(loadpath, rel_name);
            eval(['s = ' rel_name ';']);
            eval(['fs = ' rel_name '.fs;']);
            s.v = s.v(:);
            
            %% Eliminate very low frequencies
            s_filt = elim_sub_cardiac(s, up);
            nan_els = isnan(s_filt.v);
            s_filt.t = s_filt.t(~nan_els);
            s_filt.v = s_filt.v(~nan_els);
            
            %% Detect QRS Waves
            pk_inds = feval(up.al.options.RDt{current_opt_no}, s_filt, fs, subj, up);
            
            %% Refine detected QRS spikes
            % The detected peaks are currently somewhere in the middle of the R-area, so need to search for the max ecg value within the tolerance either side of the detected peak:
            max_HR_for_qrs = 300;
            tolerance = round(fs/(max_HR_for_qrs/60));
            if ~isempty(pk_inds)
                % refine the peaks and store in a new variable (ECG_PKS_inds)
                ref_pk_inds = nan(length(pk_inds),1);
                for peak_ind = 1 : length(pk_inds)                             % for each peak
                    % Find tolerance limits, centred on the current peak index:
                    if pk_inds(peak_ind) <= tolerance                          % if the peak is right at the start of the recording, set the lower tolerance limit to the start of the recording
                        lower_lim = 1;
                    else
                        lower_lim = pk_inds(peak_ind)-tolerance;
                    end
                    if (length(s.v)-pk_inds(peak_ind)) <= tolerance % if the peak is right at the end of the recording, set the upper tolerance limit to the end of the recording
                        upper_lim = length(s.v);
                    else
                        upper_lim = pk_inds(peak_ind)+tolerance;
                    end
                    % Find the maximum ecg value within the tolerance limits:
                    [~, max_ind] = max(s.v(lower_lim : upper_lim));
                    % Store the index of this maximum value, referenced to the section_data inds:
                    ref_pk_inds(peak_ind) = lower_lim-1+max_ind;
                    clear max_ind
                end
                clear peak_ind
            else
                % if no peaks were detected then give empty results:
                ref_pk_inds = [];
            end
            
            %% Eliminate any peaks which are the same
            ref_pk_inds = unique(ref_pk_inds);
            
            %% Find troughs
            % search 0.1s either side
            troughs.i = nan(length(ref_pk_inds)-1,1);
            search_min = ref_pk_inds - ceil(0.1*fs);
            search_min(search_min<1) = 1;
            search_max = ref_pk_inds + ceil(0.1*fs);
            search_max(search_max>length(s.v)) = 1;
            for peak_ind = 1 : (length(ref_pk_inds)-1)
                [~, rel_el] = min(s.v(search_min(peak_ind):search_max(peak_ind)));
                rel_el = search_min(peak_ind)-1+rel_el;
                troughs.i(peak_ind) = rel_el;
            end
            
            %% Save processed data
            eval([save_name '.fs = fs;']);
            eval([save_name '.p.v = s.v(ref_pk_inds);']);
            eval([save_name '.p.t = s.t(ref_pk_inds);']);
            eval([save_name '.tr.v = s.v(troughs.i);']);
            eval([save_name '.tr.t = s.t(troughs.i);']);
            % Identify start and end times of raw data (for resampling)
            eval([save_name '.timings.t_start = s.t(1);']);
            eval([save_name '.timings.t_end = s.t(end);']);
            save_or_append_data
            
        end
        
    end
    
end

end

function pk_inds = GC(s, fs, subj, up)

%% This uses Prof Gari Clifford's "rpeakdetect.m" script
% This script can be downloaded from:
%    http://www.mit.edu/~gari/CODE/ECGtools/ecgBag/rpeakdetect.m
%
% The following is an excerpt from the script:
%
% Written by G. Clifford gari@ieee.org and made available under the 
% GNU general public license. If you have not received a copy of this 
% license, please download a copy from http://www.gnu.org/
%
% Please distribute (and modify) freely, commenting
% where you have added modifications. 
% The author would appreciate correspondence regarding
% corrections, modifications, improvements etc.

[~, ~, ~, pk_inds, ~, ~]  = rpeakdetect(s.v(:),fs);

end