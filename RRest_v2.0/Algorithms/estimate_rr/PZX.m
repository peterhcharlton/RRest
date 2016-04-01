function rr = PZX(rel_data, wins, up)
%PZX estimates the RR using a combination of zero-crossing detection, and
% peak detection.
%	            PZX(option, up)
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

%% Cycle through windows
rr.t = mean([wins.t_start(:)' ; wins.t_end(:)']); rr.t = rr.t(:);
rr.v = nan(length(rr.t),1);

for win_no = 1 : length(wins.t_start)
    
    % extract relevant data
    rel_els = find(rel_data.t >= wins.t_start(win_no) & rel_data.t < wins.t_end(win_no));
    data.v = rel_data.v(rel_els);
    data.t = rel_data.t(rel_els);
    
    good_els = ~isnan(data.v);
    data.v = data.v(good_els); data.v = data.v(:);
    data.t = data.t(good_els); data.t = data.t(:);
    
    %% identify individual breaths from the raw signal using 3-pt peak detection:
    diffs_on_left_of_pt = diff(data.v); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt>0);
    diffs_on_right_of_pt = diff(data.v); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt<0);
    peaks = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
    
    %% identify individual breaths from the raw signal using 3-pt trough detection:
    diffs_on_left_of_pt = diff(-1*data.v); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt>0);
    diffs_on_right_of_pt = diff(-1*data.v); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt<0);
    troughs = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
    
    %% eliminate peaks or troughs which aren't in their half (above or below mean)
    mean_val = mean(data.v);
    temp_log = data.v(peaks) > mean_val;
    peaks = peaks(temp_log);
    temp_log = data.v(troughs) < mean_val;
    troughs = troughs(temp_log);
    
    %% Eliminate peaks or troughs which occur within 0.5s of the previous peak or trough
    thresh = 0.5;
    % peaks
    diffs = diff(data.t(peaks));
    rel_diffs = find(diffs > thresh); rel_diffs = rel_diffs(:);
    if isempty(rel_diffs)
        rr.v(win_no) = nan;
        continue
    end
    rel_diffs = [1; rel_diffs+1];
    peaks = peaks(rel_diffs);
    % troughs
    diffs = diff(data.t(troughs));
    rel_diffs = find(diffs > thresh); rel_diffs = rel_diffs(:);
    if isempty(rel_diffs)
        rr.v(win_no) = nan;
        continue
    end
    rel_diffs = [1; rel_diffs+1];
    troughs = troughs(rel_diffs);
    
    %% Eliminate peaks or troughs which aren't followed by their opposite
    comp.t = [peaks; troughs];
    comp.log = [ones(length(peaks),1); zeros(length(troughs),1)];
    [comp.t, orders] = sort(comp.t);
    comp.log = comp.log(orders);
    diffs = diff(comp.log); diffs = [1; diffs];
    comp.t = comp.t(diffs ~= 0);
    comp.log = comp.log(diffs ~= 0);
    peaks = comp.t(comp.log==1);
    troughs = comp.t(comp.log==0);    
    
    %% Calc RR
    breaths = peaks;
    % Only using time period spanning a whole number of breaths
    win_length = data.t(breaths(end)) - data.t(breaths(1));
    rr.v(win_no) = 60*(length(breaths)-1)/win_length;
    
end

end