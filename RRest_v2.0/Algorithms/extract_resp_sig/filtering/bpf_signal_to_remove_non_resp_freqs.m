function s_filt = bpf_signal_to_remove_non_resp_freqs(s, Fs, up)
% Filter pre-processed signal to remove freqs outside of the range of plausible resp freqs

s.fs = Fs;

%% Window signal to reduce edge effects
s_win = tukey_win_data(s,up);

%% HPF to eliminate VLFs
s_evlfs = s_win;
s_evlfs.v = elim_vlfs(s_win, up);

%% LPF to remove freqs above resp
s_filt.t = s_evlfs.t;
s_filt.v = lp_filter_signal_to_remove_freqs_above_resp(s_evlfs.v, s_evlfs.fs, up);
s_filt.fs = s_evlfs.fs;

end