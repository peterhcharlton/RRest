function calculate_resp_sig_cc(up)
%CALCULATE_RESP_SIG_CC calculates correlation coefficients (CCs) between the
% extracted respiratory signals and reference respiratory signal
%	            calculate_resp_sig_cc(up)
%
%	Inputs:
%		data            data files should be stored in the specified format
%       up              universal parameters structure
%
%	Outputs:
%       	for each subject, n:
%       n_cc.m          - a file of CC values
%

fprintf('\n--- Calculating CCs ');
%% Extract list of resp signals from the first pt:
subj = up.paramSet.subj_list(1);
loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.respSigs];
if exist(loadpath, 'file')
    filecontents = whos('-file', loadpath);
    respSigs = extractfield(filecontents, 'name');
else
    warning(['No respiratory signals found for Subject ', num2str(subj) '.'])
end

% For each subject
for subj = up.paramSet.subj_list
    
    loaded_this_subj_respSigs = 0;
    %% Make window timings if necessary
    identify_subj_wins(subj, up);
    %% Cycle through each resp signal
    for respSig_no = 1:length(respSigs)
        % Skip if this processing has been done previously
        save_name = respSigs{respSig_no};
        savepath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.cc, '.mat'];
        exist_log = check_exists(savepath, save_name);
        if exist_log
            continue
        end
        
        % load this subject's resp sig data if it hasn't yet been loaded
        if ~loaded_this_subj_respSigs
            % Signals
            load([up.paths.data_save_folder, num2str(subj), up.paths.filenames.respSigs]);
            % Window timings
            load([up.paths.data_save_folder, num2str(subj), up.paths.filenames.win_timings, '.mat']);
            % ref RR
            load([up.paths.data_save_folder, num2str(subj), up.paths.filenames.ref_rrs]);
            loaded_this_subj_respSigs = 1;
        end
        
        % load reference resp sig data if it hasn't been loaded
        if ~exist('data', 'var')
            load([up.paths.data_load_folder, up.paths.data_load_filename]);
        end
                
        %up.paramSet.ref_method = 'paw';
        
        % Extract respiratory signals
        ext_sig = extract_resp_sigs(data, subj, up);
        % Determine SNRs
        snr_sig = calc_snrs(ext_sig, subj, wins, up);
        clear ext_sig
        % LPF to exclude freqs above resp
        rel_data.ref = lpf_to_exclude_resp(snr_sig, subj, up);
        clear snr_sig
        
        % Identify the relevant respSig data
        eval(['rel_data.ext = ' respSigs{respSig_no} ';']);
        % Identify the relevant ref RR data
        rel_data.rr = rr_ref;
        
        %% Calculate CC between this resp sig and the reference respiratory signal
        if (length(rel_data.ext.t) == 1 && isnan(rel_data.ext.t)) || sum(isnan(rel_data.ext.v))==length(rel_data.ext.v)
            temp_cc.t = mean([wins.t_start(:)' ; wins.t_end(:)']); temp_cc.t = temp_cc.t(:);
            [temp_cc.MSC, temp_cc.SNR, temp_cc.CCp] = deal(nan(length(temp_cc.t),1));
        else
            temp_cc = calc_cc(rel_data, wins, up);
        end
        % store this series of CCs:
        eval([respSigs{respSig_no} ' = temp_cc;']);
        clear temp_cc
        %% Save RRs to file
        save_or_append_data
    end
    clear ekg* ppg* wins            % clear resp sigs
end

end

function cc = calc_cc(data, wins, up)

%% Setup
downsample_freq = up.paramSet.fft_resample_freq;    % it would be worth changing this - it changes the answer (try 1 Hz e.g.)

%% Sample both signals at the same points at 5 Hz

% Identify latest start time and earliest end time
rel_ref_start_el = find(data.ref.t>= wins.t_start(1), 1);
rel_ext_start_el = find(data.ext.t>= wins.t_start(1), 1);
rel_ref_end_el = find(data.ref.t<= wins.t_end(end), 1, 'last');
rel_ext_end_el = find(data.ext.t<= wins.t_end(end), 1, 'last');
t_start = max([data.ref.t(rel_ref_start_el), data.ext.t(rel_ext_start_el)]);
t_end = min([data.ref.t(rel_ref_end_el), data.ext.t(rel_ext_end_el)]);
% create re-sample grid
resample_grid = t_start:(1/downsample_freq):t_end;

for sig_type = {'ref', 'ext'}
    
    eval(['temp_data = data.' sig_type{1,1} ';']);
    
    % Downsample onto the re-sample grid
    temp_data.filt.t = resample_grid;
    temp_data.filt.v = interp1(temp_data.t, temp_data.v, temp_data.filt.t);
    temp_data.filt.v = temp_data.filt.v(:);
    temp_data.filt.t = temp_data.filt.t(:);
    
    % Remove non-resp freqs
    temp_data.filt = bpf_signal_to_remove_non_resp_freqs(temp_data.filt, downsample_freq, up);
    
    %temp_data.filt.t = downsample(temp_data.t, temp_data.fs/downsample_freq);
    %temp_data.filt.v = decimate(temp_data.v, temp_data.fs/downsample_freq);
    %temp_data.filt.v = detrend(temp_data.filt.v);
    
    eval(['data.' sig_type{1,1} '= temp_data.filt;']);
end
clear temp_data

%% Cycle through windows
rr.t = mean([wins.t_start(:)' ; wins.t_end(:)']); rr.t = rr.t(:);
rr.v = nan(length(rr.t),1);

[CCp, MSC, SNR] = deal(nan(length(wins.t_start),1));
for win_no = 1 : length(wins.t_start)
    
    % find true rr
    true_rr = data.rr.v(win_no);
    if data.rr.snr_log(win_no) == 0 || isnan(true_rr)
        continue
    end
    
    % extract relevant data for this window
    
    for sig_type = {'ref', 'ext'}
        
        eval(['curr_data = data.' sig_type{1,1} ';']);
        
        temp_data.fs = curr_data.fs;
        
        rel_els = find(curr_data.t >= wins.t_start(win_no) & curr_data.t < wins.t_end(win_no));
        temp_data.v = curr_data.v(rel_els);
        temp_data.t = curr_data.t(rel_els); clear rel_els
        
        good_els = ~isnan(temp_data.v);
        temp_data.v = temp_data.v(good_els);
        temp_data.t = temp_data.t(good_els); clear good_els
        
        eval(['rel_data.' sig_type{1,1} '= temp_data;']);
        
    end
    
    % Find the lag between the two:  
    [r,lags] = xcorr(rel_data.ref.v, rel_data.ext.v);
    % only use the lags corresponding to a max of half the maximum breath duration either side:
    max_breath_dur = 60/true_rr; % in secs
    max_no_samps = floor(max_breath_dur*downsample_freq);
    mid_lag = (length(lags)+1)/2;
    rel_lags = mid_lag - max_no_samps : mid_lag + max_no_samps;
    r = r(rel_lags);
    lags = lags(rel_lags);
    % find the relevant lag:
    [~, rel_el] = max(abs(r));
    rel_lag = lags(rel_el); % a negative lag means the second signal has peaks after the first signal, so you need to take elements off the start of the second signal, and add elements to the start of the first signal.
    
    % offset the signals accordingly:
    old = rel_data;
    if rel_lag <0
        rel_data.ref.v = [zeros(abs(rel_lag),1); rel_data.ref.v];
        rel_ref_els = [false(abs(rel_lag),1); true(length(rel_data.ref.v),1)];
        rel_data.ext.v = [rel_data.ext.v; zeros(abs(rel_lag),1)];
        rel_ext_els = [true(length(rel_data.ext.v),1); false(abs(rel_lag),1)];
    elseif rel_lag > 0
        rel_data.ext.v = [zeros(rel_lag,1); rel_data.ext.v];
        rel_ext_els = [false(abs(rel_lag),1); true(length(rel_data.ext.v),1)];
        rel_data.ref.v = [rel_data.ref.v; zeros(rel_lag,1)];
        rel_ref_els = [true(length(rel_data.ref.v),1); false(abs(rel_lag),1)];
    else
        rel_ref_els = true(length(rel_data.ref.v),1);
        rel_ext_els = true(length(rel_data.ext.v),1);
    end
    rel_els = rel_ref_els & rel_ext_els;
    rel_data.ext.v = rel_data.ext.v(rel_els);
    rel_data.ref.v = rel_data.ref.v(rel_els);
    
    % orientate the signals accordingly:
    rel_mag = r(rel_el);
    if rel_mag < 0
        rel_data.ext.v = -1*rel_data.ext.v;
    end
    
    % Find pearson's correlation coefficient:
    CCp(win_no) = corr(rel_data.ref.v, rel_data.ext.v);
    
    % Find magnitude squared coherence at the true respiratory freq:
    [MSC_spect,F] = mscohere(rel_data.ref.v, rel_data.ext.v, [], [], [], downsample_freq);
    F = F*60; % in bpm
    [~, rel_el] = min(abs(F-true_rr)); 
    MSC(win_no) = MSC_spect(rel_el);
    
%     % Find SNR:
%     SNR(win_no) = snr(old.ext.v, downsample_freq, 6);
    
    % Find Power spectral analysis: (see http://doi.org/10.1111/j.1399-6576.2007.01375.x)
    
    %  - Find the Welch periodogram of ref signal
    segLen = 2^nextpow2(12*downsample_freq);
    noverlap = segLen/2;
    [tempr.power, tempr.freqs] = pwelch(old.ref.v,segLen,noverlap, [], downsample_freq);
    %  - Find symmetric area containing 75% of power in spectrum
    [~, centre_el] = max(tempr.power);
    curr_prop_power = 0; curr_width = 1; % width in samples
    while curr_prop_power < 0.75
        rel_els = centre_el-curr_width: centre_el+curr_width;
        rel_els = rel_els(rel_els>0);
        curr_prop_power = sum(tempr.power(rel_els))/sum(tempr.power);
        curr_width = curr_width + 1;
    end
    lower_lim = tempr.freqs(rel_els(1));
    upper_lim = tempr.freqs(rel_els(end));
    % - Find the Welch periodogram of the ext signal
    [tempe.power, tempe.freqs] = pwelch(old.ext.v,segLen,noverlap, [], downsample_freq);
    rel_els = tempe.freqs>= lower_lim & tempe.freqs <= upper_lim;
    prop_ext_power = sum(tempe.power(rel_els))/sum(tempe.power);
    
    SNR(win_no) = prop_ext_power;
    
end

t = mean([wins.t_start(:)' ; wins.t_end(:)']); t = t(:);

cc.t = t;
cc.CCp = CCp;
cc.MSC = MSC;
cc.SNR = SNR;

end

function new_sig = extract_resp_sigs(data, subj, up)

%% See which resp signals are present:
if sum(strcmp(up.paramSet.ref_method(1:3), 'imp'))
    % select impedance:
    if strcmp(fieldnames(data(subj).ref.resp_sig), 'imp')
        imp.v = data(subj).ref.resp_sig.imp.v;
        imp.fs = data(subj).ref.resp_sig.imp.fs;
        imp.t = (1/imp.fs)*(1:length(imp.v));
        new_sig = imp;
    elseif strcmp(fieldnames(data(subj).ref.resp_sig), 'unknown')
        unk.v = data(subj).ref.resp_sig.unknown.v;
        unk.fs = data(subj).ref.resp_sig.unknown.fs;
        unk.t = (1/unk.fs)*(1:length(unk.v));
        new_sig = unk;
    end
elseif sum(strcmp(up.paramSet.ref_method, 'paw'))
    % select paw:
    paw.v = data(subj).ref.resp_sig.paw.v;
    paw.fs = data(subj).ref.resp_sig.paw.fs;
    paw.t = (1/paw.fs)*(1:length(paw.v));
    % store
    new_sig = paw;
end
if sum(strcmp(up.paramSet.ref_method, 'co2'))
    % select co2:
    co2.v = data(subj).ref.resp_sig.co2.v;
    co2.fs = data(subj).ref.resp_sig.co2.fs;
    co2.t = (1/co2.fs)*(1:length(co2.v));
    % store
    new_sig = co2;
elseif sum(strcmp(up.paramSet.ref_method(1:3), 'ban'))
    % select chest band:
    band.v = data(subj).ref.resp_sig.band.v;
    band.fs = data(subj).ref.resp_sig.band.fs;
    band.t = (1/band.fs)*(1:length(band.v));
    new_sig = band; 
end

end

function sig = calc_snrs(sig, subj, wins, up)

sig.snr = nan(length(wins.t_start),1);
for win_no = 1 : length(wins.t_start)
    % select data for this window
    rel_els = find(sig.t >= wins.t_start(win_no) & sig.t <= wins.t_end(win_no));
    rel_data.t = sig.t(rel_els);
    rel_data.v = sig.v(rel_els);
    
    if strcmp(up.paramSet.ref_method, 'paw')
        % calc snr
        sig.snr(win_no) = snr(rel_data.v);
    end
end
end

function lpf_sig = lpf_to_exclude_resp(sig, subj, up)

%% Window signal to reduce edge effects
duration_of_signal = sig.t(end) - sig.t(1);
prop_of_win_in_outer_regions = 2*up.paramSet.tukey_win_duration_taper/duration_of_signal;
tukey_win = tukeywin(length(sig.v), prop_of_win_in_outer_regions);
d_s_win = sig;    % copy time and fs
d_s_win.v = detrend(sig.v(:)).*tukey_win(:);

%% LPF to remove freqs above resp
lpf_sig.t = d_s_win.t;
lpf_sig.v = lp_filter_signal_to_remove_freqs_above_resp(d_s_win.v, d_s_win.fs, up, 'ref_rr');
lpf_sig.fs = d_s_win.fs;
lpf_sig.snr = sig.snr;
% NB: if you use the e_vlf signal then the freqs below resp have already been removed.

end
