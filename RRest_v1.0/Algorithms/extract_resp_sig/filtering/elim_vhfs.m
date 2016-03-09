function s_filt = elim_vhfs(s, fs, filt_characteristics)
%% Filter signal to remove VHFs

%% Eliminate nans
s(isnan(s)) = mean(s(~isnan(s)));

%% Check to see if sampling freq is at least twice the freq of interest
if (filt_characteristics.Fpass/(fs/2)) >= 1
    % then the fs is too low to perform this filtering
    s_filt = s;
    return
end

%% Create filter
% parameters for the low-pass filter to be used
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([filt_characteristics.Fstop filt_characteristics.Fpass]/(fs/2), [1 0], [filt_characteristics.Dstop filt_characteristics.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% Gives a -3 dB cutoff at cutoff_freq Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.139;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(fs/2);

%% Remove VHFs
s_dt=detrend(s);
s_filt = filtfilt(AMfilter.numerator, 1, s_dt);
end