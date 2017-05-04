function rr = FTS(rel_data, wins, up)
%FTS calculates the frequency spectrum of a signal using the FFT, and finds
% the RR.
%	            FTS(rel_data, wins, up)
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
    
    %% Now that we have a processed PPG-waveform, i.e. the 'respiratory signal', lets find the FFT
    
    % Find FFT
    WINLENGTH = length(data.filt.v);
    NFFT = 2^nextpow2(WINLENGTH);
    HAMMWIN = hamming(WINLENGTH);
    HAMMWIN = HAMMWIN(:);
    f_nyq = downsample_freq/2;
    FREQS = f_nyq.*linspace(0, 1, NFFT/2+1);            % Array of correspondent FFT bin frequencies, in BR (RPM)
    WINDATA = detrend(data.filt.v);                      % Remove the LSE straight line from the data
    WINDATA = WINDATA .* HAMMWIN;
    myFFT = fft(WINDATA, NFFT);
    myFFT = myFFT(1 : NFFT/2 + 1);
    myFFT = 2.*abs(myFFT/NFFT);
    psdx = (1/(downsample_freq*NFFT)) * abs(myFFT).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    clear myFFT WINDATA f_nyq HAMMWIN NFFT WINLENGTH
    
    data.power = 10*log10(psdx); data.power = data.power(:);
    data.freqs = FREQS; data.freqs = data.freqs(:);
    
    % Find spectral peak
    [rr.v(win_no), rr.f{win_no}, rr.p{win_no}] = find_spectral_peak(data, up);
    
    clear data
    
end

end