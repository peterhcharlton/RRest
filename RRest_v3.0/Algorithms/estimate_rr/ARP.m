function rr = ARP(rel_data, wins, up)
%ARP estimates the RR using auto-regressive all-pole modelling.
%
%	            ARP(option, up)
%
%       Based on Marco Pimentel's "armodel.m"
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
%

%% Setup
downsample_freq = up.paramSet.ar_resample_freq;    % it would be worth changing this - it changes the answer (try 1 Hz e.g.)
true_fs = rel_data.fs;

%% Cycle through windows
rr.t = mean([wins.t_start(:)' ; wins.t_end(:)']); rr.t = rr.t(:);
[rr.v, rr.pole_mag, rr.pole_angle] = deal(nan(length(rr.t),1));

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
    [a,~] = arburg(data.filt.v,up.paramSet.ar_model_order);
    if length(unique(data.filt.v)) > 1 % there is a waveform, and not a flatline
        poles=roots(a);
    else
        rr.v(win_no) = nan;
        rr.pole_mag(win_no) = 0;
        rr.pole_angle(win_no) = nan;
        continue
    end
    
    angles=angle(poles).*(downsample_freq/(2*pi))*60;   % in bpm
    % find rel poles
    rel_pole_els = angles > up.paramSet.rr_range(1) & angles < up.paramSet.rr_range(2);
    if sum(rel_pole_els) == 0
        rr.v(win_no) = nan;
        rr.pole_mag(win_no) = 0;
        rr.pole_angle(win_no) = nan;
        continue
    end
    model_details.angles=angles(rel_pole_els)';
    model_details.mag=abs(poles(rel_pole_els));
    % select the pole with the greatest magnitude as the respiratory pole
    [~, resp_pole_el] = max(model_details.mag);
    rr.v(win_no) = model_details.angles(resp_pole_el);
    rr.pole_mag(win_no) = model_details.mag(resp_pole_el);
    rr.pole_angle(win_no) = model_details.angles(resp_pole_el)./((downsample_freq/(2*pi))*60);
    
end

end