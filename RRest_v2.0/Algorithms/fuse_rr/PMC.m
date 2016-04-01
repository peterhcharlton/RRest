function rr = PMC(data, up)
%PMC fuses RR estimates using the pole magnitude criterion.

%% Find times at which all three RRs were estimated
mutual_times = intersect(intersect(data.bw.t, data.fm.t), data.am.t);
quit_log = 0;
for mod = {'bw', 'am', 'fm'}
    eval(['temp_data = data.' mod{1,1} ';']);
    [rel_times,rel_els,~] = intersect(temp_data.t,mutual_times);
    if ~sum(strcmp(fieldnames(temp_data), 'pole_mag'))
        quit_log = 1;
    else
        eval([mod{1,1} '.est = data.' mod{1,1} '.v(rel_els);']);
        eval([mod{1,1} '.mag = data.' mod{1,1} '.pole_mag(rel_els);']);
        % eval([mod{1,1} '.ang = data.' mod{1,1} '.pole_angle(rel_els);']);
    end
end
rr.v = nan(length(rel_times),1);
rr.t = rel_times;

% quit if one of the modulations wasn't present
if quit_log
    return
end

%% Find PMC
for win_no = 1 : length(mutual_times)
    [rel_ests, rel_mags] = deal([]);
    for mod = {'bw', 'am', 'fm'}
        eval(['rel_ests = [rel_ests, ' mod{1,1} '.est(win_no)];']);
        eval(['rel_mags = [rel_mags, ' mod{1,1} '.mag(win_no)];']);
    end
    %% Find RR
    [~,rel_mod_no] = max(rel_mags);
    rr.v(win_no) = rel_ests(rel_mod_no);
end

end