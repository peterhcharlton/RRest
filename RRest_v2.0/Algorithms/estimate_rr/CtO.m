function rr = CtO(rel_data, wins, up)
%CtO estimates the RR using an implementation of Count Orig, as specified in:
% A. Schäfer et al., “Estimation of breathing rate from respiratory sinus arrhythmia: comparison of various methods,” Ann. Biomed. Eng., vol. 36, no. 3, pp. 476–85, Mar. 2008.
%	            CtO(option, up)
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
    q3 = quantile(data.v(peaks), 0.75);
    thresh = 0.2*q3;
    
    %% find relevant peaks and troughs
    extrema = sort([peaks(:); troughs(:)]);
    rel_peaks = peaks(data.v(peaks) > thresh);
    rel_troughs = troughs(data.v(troughs) < 0);
    
    %% find valid breathing cycles
    % valid cycles start with a peak:
    valid_cycles = zeros(length(rel_peaks)-1,1);
    cycle_durations = nan(length(rel_peaks)-1,1);
    for peak_no = 1 : (length(rel_peaks)-1)
        
        % valid if there is only one rel_trough between this peak and the
        % next
        cycle_rel_troughs = rel_troughs(rel_troughs > rel_peaks(peak_no) & rel_troughs < rel_peaks(peak_no+1));
        if length(cycle_rel_troughs) == 1
            valid_cycles(peak_no) = 1;
            cycle_durations(peak_no) = data.t(rel_peaks(peak_no+1)) - data.t(rel_peaks(peak_no));
        end
        
    end
    
    if sum(valid_cycles) == 0 
        rr.v(win_no) = nan;
        continue
    end
    
    %% Calc RR

    % Using average breath length
    ave_breath_duration = nanmean(cycle_durations);
    rr.v(win_no) = 60/ave_breath_duration;
    
end

end