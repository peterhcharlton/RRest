function FPt(option, up)
%FPT detects Fiducial Points from PPG peak and trough annotations
% as specified in PC's literature review.
%	            FPt(option, up)
%
%	Inputs:
%		option          the option which has led to this function being used
%       up              universal parameters structure
%
%	Outputs:
%       ...
%

fprintf('\n--- Identifying Fiducial Points ');
log_int_respSig = 1;             % Has value 1 unless this is a final respiratory signal

for subj = up.paramSet.subj_list
    
    sig_type = option(1:3);
    %% Cycle through each signal of this type
    eval(['sigs = up.paramSet.' sig_type '_sigs;']);
    for sig_no = 1 : length(sigs)
        curr_sig = sigs{sig_no};
        curr_sig_type = curr_sig(1:3);
        
        %% Skip if this processing has been done previously
        iden_resp_sig_file_ending
        savepath = [up.paths.data_save_folder, num2str(subj), ending];
        filecontents = whos('-file', savepath);
        var_names = extractfield(filecontents, 'name');
        rel_log = zeros(size(var_names));
        if strcmp(sig_type, 'ekg')
            rel_string = [curr_sig, up.paths.filenames.qrss];
        else
            rel_string = [curr_sig, up.paths.filenames.pulse_peaks];
        end
        for s = 1 : length(var_names)
            if strfind(var_names{s}, rel_string)
                rel_log(s) = 1;
            end
        end
        rel_var_names = var_names(logical(rel_log));
        for rel_var_name_no = 1 : length(rel_var_names)
            temp = strfind(rel_var_names{rel_var_name_no},'_');
            start_el = temp(1); clear temp
            eval(['save_name = ''' curr_sig, up.paths.filenames.fid_pts, rel_var_names{rel_var_name_no}(start_el:end) ''';']); clear start_el
            exist_log = check_exists(savepath, save_name);
            if exist_log
                continue
            end
            
            %% Load relevant data
            loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.int_respSigs];
            
            if strcmp(curr_sig_type, 'ppg')
                rel_name1 = rel_var_names{rel_var_name_no};
                rel_name2 = [curr_sig up.paths.filenames.elim_vhf];
                load(loadpath, rel_name1, rel_name2);
                eval(['rel_data.beats = ' rel_name1 ';']);
                eval(['rel_data.s = ' rel_name2 ';']);
                eval(['rel_data.fs = ' rel_name1 '.fs;']);
                eval(['rel_data.timings = ' rel_name1 '.timings;']);
            else
                %% ECG
                rel_name1 = rel_var_names{rel_var_name_no};
                rel_name2 = [curr_sig up.paths.filenames.elim_vhf];
                load(loadpath, rel_name1, rel_name2);
                eval(['rel_data.beats = ' rel_name1 ';']);
                eval(['rel_data.s = ' rel_name2 ';']);
                eval(['rel_data.fs = ' rel_name1 '.fs;']);
                eval(['rel_data.timings = ' rel_name1 '.timings;']);
            end
            
            %% Identify Fiducial Points
            
            if strcmp(curr_sig_type, 'ppg')
                %% PPG Peaks
                % Find peaks as max between detected onsets
                temp.p_max.t = nan(length(rel_data.beats.tr.t)-1,1);
                temp.p_max.v = nan(length(rel_data.beats.tr.t)-1,1);
                for beat_no = 1 : (length(rel_data.beats.tr.t)-1)
                    rel_range = find(rel_data.s.t >= rel_data.beats.tr.t(beat_no) & rel_data.s.t < rel_data.beats.tr.t(beat_no+1));
                    [~, rel_el] = max(rel_data.s.v(rel_range));
                    temp.p_max.t(beat_no) = rel_data.s.t(rel_range(rel_el));
                    temp.p_max.v(beat_no) = rel_data.s.v(rel_range(rel_el));
                end
                
                %% PPG Troughs
                
                % Find troughs as min between detected peaks
                temp.tr_min.t = nan(length(temp.p_max.t)-1,1);
                temp.tr_min.v = nan(length(temp.p_max.t)-1,1);
                for beat_no = 1 : (length(temp.p_max.t)-1)
                    rel_range = find(rel_data.s.t >= temp.p_max.t(beat_no) & rel_data.s.t < temp.p_max.t(beat_no+1));
                    [~, rel_el] = min(rel_data.s.v(rel_range));
                    temp.tr_min.t(beat_no) = rel_data.s.t(rel_range(rel_el));
                    temp.tr_min.v(beat_no) = rel_data.s.v(rel_range(rel_el));
                end
                
            else
                %% ECG Peaks
                % Used to find peaks as max between detected onsets, but it now
                % appears that this doesn't work if it detects the onsets as a
                % mixture of Q and S waves.
                
                temp.p_max = rel_data.beats.p;
                
                %% ECG Troughs
                % Find troughs as min within search range 0.1s before peaks
                temp.tr_min.t = nan(length(temp.p_max.t),1);
                temp.tr_min.v = nan(length(temp.p_max.t),1);
                thresh = 0.1;
                for beat_no = 1 : (length(temp.p_max.t))
                    rel_range = find(rel_data.s.t >= (temp.p_max.t(beat_no) - thresh) & rel_data.s.t < temp.p_max.t(beat_no));
                    % used to be:     rel_range = find(rel_data.s.t >= (temp.p_max.t(beat_no) - thresh) & rel_data.s.t < (temp.p_max.t(beat_no) + thresh));
                    [~, rel_el] = min(rel_data.s.v(rel_range));
                    if ~isempty(rel_el)   % it is empty if the peak is at the very first element of the signal
                        temp.tr_min.t(beat_no) = rel_data.s.t(rel_range(rel_el));
                        temp.tr_min.v(beat_no) = rel_data.s.v(rel_range(rel_el));
                    end
                end
                % get rid of any nans (arise if the peak is at the very first element of the signal)
                bad_els = isnan(temp.tr_min.t);
                temp.tr_min.t = temp.tr_min.t(~bad_els);
                temp.tr_min.v = temp.tr_min.v(~bad_els);
                
                % very ocassionally it picks out the same trough or peak twice (if two consecutive peaks are ridiculously close together so share some of the same search range)
                [~, rel_els, ~] = unique(temp.tr_min.t);
                temp.tr_min.t = temp.tr_min.t(rel_els);
                temp.tr_min.v = temp.tr_min.v(rel_els);
                
            end
            
            % Carry forward detected peaks and onsets
            temp.det_p.t = rel_data.beats.p.t;
            temp.det_p.v = rel_data.beats.p.v;
            temp.det_tr.t = rel_data.beats.tr.t;
            temp.det_tr.v = rel_data.beats.tr.v;
            
            % carry forward fs and timings
            temp.fs = rel_data.fs;
            temp.timings = rel_data.timings;
            
            %% Save processed data
            eval([save_name ' = temp;']);
            save_or_append_data
            
        end
        
    end
    
end

end