%% iden_resp_sig_file_ending - finds out what the appropriate file ending is for this resp sig file
% log_int_respSig - a logical as to whether (1) this is an intermediate
% resp sig file or not(0).

if log_int_respSig
    ending = up.paths.filenames.int_respSigs;
else
    ending = up.paths.filenames.respSigs;
end