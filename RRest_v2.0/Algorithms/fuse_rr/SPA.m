function rr = SPA(data, up)
%SPA fuses RR estimates using Spectral Peak Averaging, modified from the method described in:
% Lázaro Plaza, J., 2015. Non-invasive techniques for respiratory information extraction based on pulse photoplethysmogram and electrocardiogram. Universidad Zaragoza.
% to remove the need for a previous reference RR
eta = 0.65;
lambda = 0.05;
delta = 0.1;

bw_est = data.bw;
am_est = data.am;
fm_est = data.fm;

%% Cycle through times at which RR was estimated
rr.v = nan(length(bw_est.v),1);
rr.t = bw_est.t;
for win_no = 1 : length(rr.v)
    
    % extract rel data
    rrEst_values = [bw_est.v(win_no), am_est.v(win_no), fm_est.v(win_no)];
    
    % eliminate any spectra which had rr ests (i.e. max peaks) outside of freq range:
    rel_els = ~isnan(rrEst_values); clear rrEst_values
    
    % extract spectra (taking care not to extract those which are nans and
    % so don't have freq spectra)
    spectPowers = [];
    spectFreqs = [];
    if rel_els(1)
        spectPowers = [spectPowers, bw_est.p{win_no}];
        spectFreqs = [spectFreqs, bw_est.f{win_no}];
    end
    if rel_els(2)
        spectPowers = [spectPowers, am_est.p{win_no}];
        spectFreqs = [spectFreqs, am_est.f{win_no}];
    end
    if rel_els(3)
        spectPowers = [spectPowers, fm_est.p{win_no}];
        spectFreqs = [spectFreqs, fm_est.f{win_no}];
    end
    
    % check freqs are the same:
    if sum(sum(spectFreqs, 2) == (size(spectFreqs,2)*spectFreqs(:,1))) ~= length(spectFreqs)
        error('check this')
    end
    
    % normalise spectra
    for spect_no = 1 : size(spectPowers,2)
        spectPowers(:,spect_no) = spectPowers(:,spect_no)./sum(spectPowers(:,spect_no));
    end
    clear spect_no
    
    % find power props for each spectrum
    p = find_power_props(spectPowers, bw_est, delta, win_no, up);
    
    % evaluate criteria
    rel_els = eval_criteria(p, eta, lambda);
    if isempty(rel_els)
        % repeat with twice delta
        p = find_power_props(spectPowers, bw_est, 2*delta, win_no, up);
        rel_els = eval_criteria(p, eta, lambda);
    end
    if isempty(rel_els)
        rr.v(win_no) = nan;
    else
        % find mean spectrum
        s_mean.power = sum(spectPowers(:,rel_els),2); s_mean.power = s_mean.power./sum(s_mean.power);
        s_mean.freqs = bw_est.f{win_no};
        
        % Find the respiratory peak
        
        freq_range = up.paramSet.rr_range/60;
        
        cand_els = zeros(length(s_mean.power),1);
        for s = 2 : (length(s_mean.power)-1)
            if s_mean.power(s) > s_mean.power(s-1) & s_mean.power(s) > s_mean.power(s+1) & s_mean.freqs(s) > freq_range(1) & s_mean.freqs(s) < freq_range(2)
                cand_els(s) = 1;
            end
        end
        clear s
        cand_els = find(cand_els);
        
        [~, r_el] = max(s_mean.power(cand_els));
        r_el = cand_els(r_el);
        r_freq = s_mean.freqs(r_el);
        if ~isempty(r_freq)
            rr.v(win_no) = 60*r_freq;
        else
            rr.v(win_no) = nan;
        end
        clear r_freq r_el cand_els freq_range s_mean
    end
    clear p
end

end

function p = find_power_props(spectPowers, bw_est, delta, win_no, up)

for spect_no = 1 : size(spectPowers,2)
    spect.p = spectPowers(:,spect_no);
    spect.f = bw_est.f{win_no};
    % find locations of max peaks
    [~, max_el] = max(spect.p);
    peak_f = spect.f(max_el);
    % find spectral power in interval centred on peak:
    lower_lim = max([up.paramSet.rr_range(1)/60, peak_f-(0.4*delta)]);
    upper_lim = min([up.paramSet.rr_range(2)/60, peak_f+(0.4*delta)]);
    total_lower_lim = peak_f-delta;
    total_upper_lim = peak_f+delta;
    p(spect_no) = sum(spect.p(spect.f >= lower_lim & spect.f <= upper_lim))/ ...
        sum(spect.p(spect.f >= total_lower_lim & spect.f <= total_upper_lim));
    clear lower_lim upper_lim total_lower_lim total_upper_lim peak_f max_el spect
end

end

function rel_els = eval_criteria(p, eta, lambda)
    chi_1 = p > eta;
    chi_2 = p > (max(p) - lambda);
    rel_els = chi_1 & chi_2;
end