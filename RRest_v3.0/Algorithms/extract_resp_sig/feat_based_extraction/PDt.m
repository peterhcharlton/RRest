function PDt(option, up)
%PDT detects PPG pulse peaks in PPG signals.
%	            PDt(option, up)
%
%	Inputs:
%		option          the option which has led to this function being used
%       up              universal parameters structure
%
%	Outputs:
%       ...
%

fprintf('\n--- Detecting Pulse Peaks ');
log_int_respSig = 1;             % Has value 1 unless this is a final respiratory signal

for subj = up.paramSet.subj_list
    
    sig_type = option(1:3);
    %% Cycle through each signal of this type
    eval(['sigs = up.paramSet.' sig_type '_sigs;']);
    for sig_no = 1 : length(sigs)
        curr_sig = sigs{sig_no};
        
        %% Skip if this processing has been done previously
        iden_resp_sig_file_ending
        savepath = [up.paths.data_save_folder, num2str(subj), ending];
        filecontents = whos('-file', savepath);
        var_names = extractfield(filecontents, 'name');
        temp = ~cellfun(@isempty, strfind(var_names, [curr_sig, up.paths.filenames.elim_vhf]));
        rel_var_names = var_names(temp);
        for rel_var_name_no = 1 : length(rel_var_names)
            for current_opt_no = 1 : length(up.al.options.PDt)
                temp = strfind(rel_var_names{rel_var_name_no},'_');
                start_el = temp(1); clear temp
                eval(['save_name = ''' curr_sig, up.paths.filenames.pulse_peaks, up.al.options.PDt{current_opt_no}, rel_var_names{rel_var_name_no}(start_el:end) ''';']); clear start_el
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
                
                %% High-pass filter data for peak detection
                s_filt = elim_sub_cardiac(s, up);
                
                %% Detect Pulse Peaks in HPF'd version
                [peaks,onsets] = feval(up.al.options.PDt{current_opt_no}, s_filt, fs, up);
                
                eval([save_name '.fs = fs;']);
                % Identify Pulse Peaks in original version
                eval([save_name '.p.v = s.v(peaks);']);
                eval([save_name '.p.t = s.t(peaks);']);
                eval([save_name '.tr.v = s.v(onsets);']);
                eval([save_name '.tr.t = s.t(onsets);']);
                % Identify start and end times of raw data (for resampling)
                eval([save_name '.timings.t_start = s.t(1);']);
                eval([save_name '.timings.t_end = s.t(end);']);
                
                %% Save processed data
                save_or_append_data
                
                clear peaks onsets s_filt s fs temp
                
            end
        end
        
    end
    
end

end

function [peaks,onsets] = IMS(s_filt,fs,up)

[peaks,onsets,artifs] = adaptPulseSegment(s_filt.v,fs);

end