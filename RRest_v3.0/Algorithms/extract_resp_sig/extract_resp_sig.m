function extract_resp_sig(up)
% EXTRACT_RESP_SIG extracts respiratory signals from ECG and PPG signals using
% each of the chosen techniques.
%
%               extract_resp_sig(up)
%
%	Inputs:
%		up      - universal parameters structure
%
%	Outputs:
%           for each subject, n:
%       n_int_respSigs.m     - a file of intermediate respiratory signals
%       n_respSigs.m         - a file of respiratory signals for analysis
%

%% Identify which options have sub components
sub_comps = fieldnames(up.al.sub_components);

%% Cycle through each option

for option_no = 1 : length(up.al.options.extract_resp_sig)
    
    % IF this option has sub components, perform each one in turn
    if sum(strcmp(sub_comps, up.al.options.extract_resp_sig{option_no}))
        eval(['rel_sub_comps = up.al.sub_components.'  up.al.options.extract_resp_sig{option_no} ';']);
        for sub_comp_no = 1 : length(rel_sub_comps)
            feval(rel_sub_comps{sub_comp_no}, up.al.options.extract_resp_sig{option_no}, up);
        end
    else
        % OTHERWISE perform this option
        feval(up.al.options.extract_resp_sig{option_no}, up);
    end
end

end