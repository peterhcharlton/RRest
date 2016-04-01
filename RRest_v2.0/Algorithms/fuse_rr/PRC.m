function rr = PRC(data, up)
%PRC fuses RR estimates using the pole ranking criterion.

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
        eval([mod{1,1} '.ang = data.' mod{1,1} '.pole_angle(rel_els);']);
    end
end

rr.v = nan(length(rel_times),1);
rr.t = rel_times;

% quit if one of the modulations wasn't present
if quit_log
    return
end

%% Find PRC
for win_no = 1 : length(mutual_times)
    [rel_ests, rel_mags, rel_angs] = deal([]);
    rel_mods = cell(0);
    for mod = {'bw', 'am', 'fm'}
        eval(['rel_ests = [rel_ests, ' mod{1,1} '.est(win_no)];']);
        eval(['rel_mags = [rel_mags, ' mod{1,1} '.mag(win_no)];']);
        eval(['rel_angs = [rel_angs, ' mod{1,1} '.ang(win_no)];']);
        rel_mods{length(rel_mods)+1,1} = mod{1,1};
    end
    %% Find RR
    comps = {'bwam', 'bwfm', 'amfm'};
    [prc.v, prc.rr] = deal(nan(length(comps),1));
    for comp_no = 1:length(comps)
        comp = comps(comp_no);
        % identify corresponding els
        temp1 = strfind(rel_mods, comp{1,1}(1:2));
        temp2 = strfind(rel_mods, comp{1,1}(3:4));
        rel_els = sort([find(~cellfun(@isempty,temp1)), find(~cellfun(@isempty,temp2))]);
        % find PRC for this pair
        prc.v(comp_no) = mean(rel_mags(rel_els))/abs(diff(rel_angs(rel_els)));
        prc.rr(comp_no) = mean(rel_ests(rel_els));
    end
    
    % identify relevant RR est
    [~,rel_comp_no] = max(prc.v);
    rr.v(win_no) = prc.rr(rel_comp_no);
end

end