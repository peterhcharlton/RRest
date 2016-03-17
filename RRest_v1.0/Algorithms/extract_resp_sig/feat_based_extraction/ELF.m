function ELF(option, up)
%ELF eliminates very low frequencies from resampled respiratory
% signals.
%	            ELF(option, up)
%
%	Inputs:
%		option          the option which has led to this function being used
%       up              universal parameters structure
%

fprintf('\n--- Eliminating VLFs ');
wave_type = option(1:3);

for subj = up.paramSet.subj_list
    %% Skip if this processing has been done previously
    log_int_respSig = 0;             % Has value 1 unless this is a final respiratory signal
    iden_resp_sig_file_ending
    savepath = [up.paths.data_save_folder, num2str(subj), ending];
    log_int_respSig = 1;             % Has value 1 unless this is a final respiratory signal
    iden_resp_sig_file_ending
    initloadpath = [up.paths.data_save_folder, num2str(subj), ending];
    filecontents = whos('-file', initloadpath);
    var_names = extractfield(filecontents, 'name');
    rel_log = zeros(size(var_names));
    for s = 1 : length(var_names)
        if strfind(var_names{s}, [option(1:3), up.paths.filenames.resampler])
            rel_log(s) = 1;
        end
    end
    rel_var_names = var_names(logical(rel_log));
    for rel_var_name_no = 1 : length(rel_var_names)
        eval(['save_name = ''' option(1:3), up.paths.filenames.elim_vlf2, rel_var_names{rel_var_name_no}(4:end) ''';']);
        exist_log = check_exists(savepath, save_name);
        if exist_log
            continue
        end
        
        %% Load relevant data
        rel_name = rel_var_names{rel_var_name_no};
        loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.int_respSigs];
        load(loadpath, rel_name);
        eval(['old_data = ' rel_name ';']);
        
        %% Eliminate VLFs
        data.hpfilt.t = old_data.t;
        try
            data.hpfilt.v = elim_vlfs(old_data, up);
        catch
            % if there aren't enough points to use the filter, simply carry forward the previous data
            data.hpfilt.v = old_data.v;
        end
        data.hpfilt.fs = old_data.fs;
        eval([save_name ' = data.hpfilt;']);
        
        %% Save processed data
        save_or_append_data
    end
end

end