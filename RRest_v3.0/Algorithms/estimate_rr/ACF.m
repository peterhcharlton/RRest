function rr = ACF(rel_data, wins, up)
%ACF calculates the autocorrelation function, and finds the RR.
%
%	            ACF(rel_data, wins, up)
%
%	Inputs:
%       rel_data    .t  vector of times
%                   .v  vector of resp Sig values
%                   .fs sampling freq
%       up              universal parameters structure
%       wins        .t  vector of start times
%
%	Outputs:
%       rr          .t  vector of times of estimated RRs
%                   .v  vector of estimated RRs
%                   .f  vector of freqs of power spectrum
%                   .p  vector of powers of power spectrum
%

%% Setup

downsample_freq = up.paramSet.fft_resample_freq;    % it would be worth changing this - it changes the answer (try 1 Hz e.g.)
true_fs = rel_data.fs;

%% Cycle through windows
rr.t = mean([wins.t_start(:)' ; wins.t_end(:)']); rr.t = rr.t(:);
rr.v = nan(length(rr.t),1);
rr.p = cell(length(wins.t_start),1);
rr.f = cell(length(wins.t_start),1);

for win_no = 1 : length(wins.t_start)
    
    % extract relevant data
    rel_els = find(rel_data.t >= wins.t_start(win_no) & rel_data.t < wins.t_end(win_no));
    data.v = rel_data.v(rel_els);
    data.t = rel_data.t(rel_els);
    
    good_els = ~isnan(data.v);
    data.v = data.v(good_els);
    data.t = data.t(good_els);
    
    % Downsample
    data.filt.t = downsample(data.t, true_fs/downsample_freq);
    data.filt.v = decimate(data.v, true_fs/downsample_freq);
    data.filt.v = detrend(data.filt.v);
    
    data.filt.v = data.filt.v(:);
    data.filt.t = data.filt.t(:);
    
    % Find the cross correlation
    
    [r,lags] = xcorr(data.filt.v);
    
    % find the fft of the cross correlation
    
    % Find FFT
    WINLENGTH = length(r);      % data.filt is a column vector, non?
    NFFT = 2^nextpow2(WINLENGTH);           % Size of the FFT
    HAMMWIN = hamming(WINLENGTH);          % Gots me my Hamming window.  Mmm... Hamm.
    HAMMWIN = HAMMWIN(:);
    f_nyq = downsample_freq/2;                           % Nyquist frequency.  "It's a Nyquist thing."
    FREQS = f_nyq.*linspace(0, 1, NFFT/2+1);            % Array of correspondent FFT bin frequencies, in BR (RPM)
    % NB: the 60 stays regardless of window length because this is used to calculate the freq in bpm
        
    WINDATA = detrend(r);                      % Remove the LSE straight line from the data
    % used to be: WINDATA = detrend(data.filt.v, 'constant');                      % Remove the LSE straight line from the data
    WINDATA = WINDATA .* HAMMWIN;                           % Apply the Hamm
    myFFT = fft(WINDATA, NFFT);                             % Congratulations, madame - it's an FFT
    myFFT = myFFT(1 : NFFT/2 + 1);                          % Discrard the upper half, which reflects the lower half
    myFFT = 2.*abs(myFFT/NFFT);                             % Single-sided amplitude spectrum
    psdx = (1/(downsample_freq*NFFT)) * abs(myFFT).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    
    data.power = 10*log10(psdx); data.power = data.power(:);
    data.freqs = FREQS; data.freqs = data.freqs(:);
    
    %% Find spectral peak
    [rr.v(win_no), rr.f{win_no}, rr.p{win_no}] = find_spectral_peak(data, up);
    
end

end