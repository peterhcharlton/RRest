function rr = ARS(rel_data, wins, up)
%ARS calculates the power spectral density of a signal using
% auto-regressive modelling, and finds the RR
%
%	            ARS(option, up)
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
downsample_freq = up.paramSet.ar_resample_freq;    % it would be worth changing this - it changes the answer (try 1 Hz e.g.)
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
    
    % AR modelling
    [Yy,Xx]=pburg(data.filt.v,up.paramSet.ar_model_order,[],downsample_freq);
    data.freqs=Xx;
    data.power=Yy;
    
    % Find spectral peak
    [rr.v(win_no), rr.f{win_no}, rr.p{win_no}] = find_spectral_peak(data, up);
    
    clear data
    
end

end