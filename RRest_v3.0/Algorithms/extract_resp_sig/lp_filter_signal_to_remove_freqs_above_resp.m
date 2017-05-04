function s_filt = lp_filter_signal_to_remove_freqs_above_resp(s, Fs, up, ref_rr)
%% Filter pre-processed signal to remove freqs above resp

% Find out whether this is processing ref RRs, because if so then the upper
% limit needs to be set to 60 Hz
if nargin == 4 && strcmp(ref_rr, 'ref_rr')
    % parameters for the low-pass filter to be used
    flag  = 'scale';
    Dpass = up.paramSet.elim_hf_ref_rr.Dpass;
    Dstop = up.paramSet.elim_hf_ref_rr.Dstop;
    Fstop = up.paramSet.elim_hf_ref_rr.Fstop;
    Fpass = up.paramSet.elim_hf_ref_rr.Fpass;
else
    % parameters for the low-pass filter to be used
    flag  = 'scale';
    Dpass = up.paramSet.elim_hf.Dpass;
    Dstop = up.paramSet.elim_hf.Dstop;
    Fstop = up.paramSet.elim_hf.Fstop;
    Fpass = up.paramSet.elim_hf.Fpass;
end

% create filter
[N,Wn,BETA,TYPE] = kaiserord([Fstop Fpass]/(Fs/2), [1 0], [Dstop Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.3998;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(Fs/2);

% Prepare signal
s_dt=detrend(s);
s_filt = filtfilt(AMfilter.numerator, 1, s_dt);
end