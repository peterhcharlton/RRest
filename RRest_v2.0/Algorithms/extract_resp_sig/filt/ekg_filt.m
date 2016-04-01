function ekg_filt(up)
%EKG_FILT extracts respiratory signals using various filtering techniques 
% from the ECG signal as specified in PC's literature review.
%	            ekg_filt(option, up)
%
%	Inputs:
%		option          the option which has led to this function being used
%       up              universal parameters structure
%
%	Outputs:
%       ... 
%

fprintf('\n--- Extracting Respiratory Signals using Filtering Techniques ');
log_int_respSig = 0;             % Has value 1 unless this is a final respiratory signal

for subj = up.paramSet.subj_list
    
    %% Cycle through each method
    for filt_no = 1 : length(up.al.options.ekg_filt)
        
        %% Skip if this processing has been done previously
        eval(['save_name = ''ekg', up.paths.filenames.filt '_' up.al.options.ekg_filt{filt_no} ''';']);
        iden_resp_sig_file_ending
        savepath = [up.paths.data_save_folder, num2str(subj), ending];
        exist_log = check_exists(savepath, save_name);
        if exist_log
            continue
        end
        
        %% Load relevant data
        if ~exist('data', 'var')
            load([up.paths.data_load_folder, up.paths.data_load_filename]);
        end
        % Extract EKG data
        rel_data = data(subj).ekg;
        rel_data.fs = data(subj).ekg.fs;
        
        %% Filter the raw signal using this method
        respWave = feval(up.al.options.ekg_filt{filt_no}, rel_data, up);
        
        %% Band-pass filter
        filtered_data = bpf_signal_to_remove_non_resp_freqs(respWave, respWave.fs, up);
        eval([save_name ' = filtered_data;']);
        
        %% Save processed data
        save_or_append_data
    end
    
end

end