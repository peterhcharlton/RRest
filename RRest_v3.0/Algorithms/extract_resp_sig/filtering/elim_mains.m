function s_filt = elim_mains(s, fs, filt_characteristics)
%% Filter signal to remove mains freqs

% Eliminate NaNs
s(isnan(s)) = mean(s(~isnan(s)));
% Design Filter
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord(filt_characteristics.fcuts/(fs/2), [0 1 0], [filt_characteristics.Dpass, filt_characteristics.Dstop, filt_characteristics.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);
% Detrend Signal
s_dt=detrend(s);
% Filter Signal
s_filt = filtfilt(AMfilter.numerator, 1, s_dt);
s_filt = s - s_filt;
end