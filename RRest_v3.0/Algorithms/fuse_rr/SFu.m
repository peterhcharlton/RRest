function rr = SFu(data, up)
%SFu fuses RR estimates using an implementation of smart fusion.

bw_est = data.bw.v;
am_est = data.am.v;
fm_est = data.fm.v;

%% Find times at which all three RRs were estimated
mutual_times = intersect(intersect(data.bw.t, data.fm.t), data.am.t);
for mod = {'bw', 'am', 'fm'}
    eval(['temp_data = data.' mod{1,1} ';']);
    [rel_times,rel_els,~] = intersect(temp_data.t,mutual_times);

    eval([mod{1,1} '_est = data.' mod{1,1} '.v(rel_els);']);
end

%% Cycle through times at which RR was estimated
% Find out min length incase one has fewer estimates.
min_length = min([length(bw_est), length(am_est), length(fm_est)]);
rr.v = nan(min_length,1);
rr.t = rel_times;
for s = 1 : min_length
    
    rrEst_values = [bw_est(s), am_est(s), fm_est(s)];
    
    rel_els = ~isnan(rrEst_values);
    
    standev = std(rrEst_values);
    
    if sum(rel_els) < 3 || standev > 4
        rr.v(s) = nan;
    else
        rr.v(s) = mean(rrEst_values);
    end
    
end

end