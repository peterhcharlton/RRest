function rr = CtA(rel_data, wins, up)
%CtA estimates the RR using an implementation of Count Adv, as specified in:
% A. Schäfer et al., “Estimation of breathing rate from respiratory sinus arrhythmia: comparison of various methods,” Ann. Biomed. Eng., vol. 36, no. 3, pp. 476–85, Mar. 2008.
%	            CtA(option, up)
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

true_fs = rel_data.fs;

%% Cycle through windows
rr.t = mean([wins.t_start(:)' ; wins.t_end(:)']); rr.t = rr.t(:);
rr.v = nan(length(rr.t),1);

for win_no = 1 : length(wins.t_start)
    
    %% extract relevant data
    rel_els = find(rel_data.t >= wins.t_start(win_no) & rel_data.t < wins.t_end(win_no));
    data.v = rel_data.v(rel_els);
    data.t = rel_data.t(rel_els);
    
    %% eliminate nans
    good_els = ~isnan(data.v);
    data.v = data.v(good_els); data.v = data.v(:);
    data.t = data.t(good_els); data.t = data.t(:);

    %% identify peaks
    diffs_on_left_of_pt = diff(data.v); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt>0);
    diffs_on_right_of_pt = diff(data.v); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt<0);
    peaks = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
    
    %% identify troughs
    diffs_on_left_of_pt = diff(data.v); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt<0);
    diffs_on_right_of_pt = diff(data.v); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt>0);
    troughs = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
    
    %% define threshold
    extrema = sort([peaks(:); troughs(:)]);
    amps = data.v(extrema);
    amp_diffs = abs(amps(2:end)-amps(1:(end-1)));
    q3 = quantile(amp_diffs, 0.75);
    thresh = 0.3*q3;
    
    eliminating_pairs_of_extrema = 1;
    while eliminating_pairs_of_extrema
        if length(extrema) < 3
            eliminating_pairs_of_extrema = 0;
            continue
        end
        amps = data.v(extrema);
        amp_diffs = abs(amps(2:end)-amps(1:(end-1)));
        [min_amp_diff, el] = min(amp_diffs);
        if min_amp_diff > thresh
            eliminating_pairs_of_extrema = 0;
        else
            pair_w_smallest_amp_diff = [el, el+1];
            extrema = [extrema(1:(pair_w_smallest_amp_diff(1)-1)); extrema((1+pair_w_smallest_amp_diff(2)):end)];
        end
    end
    
    if length(extrema) < 3 
        rr.v(win_no) = nan;
        continue
    end
    
    %% truncate to start and end at a maximum
    if data.v(extrema(1)) < data.v(extrema(2))
        extrema = extrema(2:end);
    end
    if data.v(extrema(end)) < data.v(extrema(end-1))
        extrema = extrema(1:end-1);
    end
    
    if length(extrema) < 3 
        rr.v(win_no) = nan;
        continue
    end
    
    %% Calc RR
    
    % find no breaths
    no_breaths = (length(extrema)-1)/2;

    % total breath duration
    total_breath_duration = data.t(extrema(end))-data.t(extrema(1));
    ave_breath_duration = total_breath_duration/no_breaths;
    rr.v(win_no) = 60/ave_breath_duration;
    
end

end