function rr = PKS(rel_data, wins, up)
%PKS estimates the RR using peak detection.
%
%	            PKS(option, up)
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
    breaths = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
    if isempty(breaths)
        rr.v(win_no) = nan;
        continue
    end
    
    %% Calc RR
    % Only using time period spanning a whole number of breaths
    win_length = data.t(breaths(end)) - data.t(breaths(1));
    rr.v(win_no) = 60*(length(breaths)-1)/win_length;
    
end

end