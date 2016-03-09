function s_filt = elim_sub_cardiac(old_data, up)
%% Filter pre-processed signal to remove frequencies below cardiac freqs

fs = old_data.fs;
s = old_data;

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Downsample
d_s = downsample_data(s, up);

%% Make filter
flag  = 'scale';        % Sampling Flag
[N,Wn,BETA,TYPE] = kaiserord([up.paramSet.elim_sub_cardiac.Fstop up.paramSet.elim_sub_cardiac.Fpass]/(d_s.fs/2), [1 0], [up.paramSet.elim_sub_cardiac.Dstop up.paramSet.elim_sub_cardiac.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);   % Calculate the coefficients using the FIR1 function.
AMfilter = dfilt.dffir(b);

%% Check frequency response
% Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.04;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(d_s.fs/2);

temp_filt = filtfilt(AMfilter.numerator, 1, d_s.v);

%% Resample
s_filt_rs.v = interp1(d_s.t, temp_filt, s.t);
s_filt.v = s.v(:)-s_filt_rs.v(:);
s_filt.t = s.t;
s_filt.fs = s.fs;
end