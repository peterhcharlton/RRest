function estimate_ref_rr(up)
%ESTIMATE_REF_RR estimates a reference RR from impedance and paw signals
%	            estimate_ref_rr(up)
%
%	Inputs:
%		data            data files should be stored in the specified format
%       up              universal parameters structure
%
%	Outputs:
%       ...
%

fprintf('\n--- Estimating Reference RRs ');

% For each subject
for subj = up.paramSet.subj_list
    
    % Skip if this processing has been done previously
    save_name = 'rr_ref';
    savepath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.ref_rrs, '.mat'];
    exist_log = check_exists(savepath, save_name);
    if exist_log
        %continue
    end
    
    %% Load window timings
    loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.win_timings, '.mat'];
    load(loadpath);
    
    %% Load data
    % if the data matrix hasn't already been loaded, load it
    if ~exist('data', 'var')
        load([up.paths.data_load_folder, up.paths.data_load_filename]);
    end
    
    %% Find additional RRs from simultaneous impedance numerics
    if up.analysis.imp_stats
        rr_ref.imp = find_rr_ref_from_ref_rrs(wins, data, subj);
    end
    
    %% Find ref RRs from breath timings
    if strcmp(up.paramSet.ref_method, 'breaths')
        rr_ref.timings = find_rr_ref_from_ref_breaths(wins, data, subj, up);
        save_or_append_data
        continue
    end
    
    %% Find ref RRs from original simultaneous RRs
    if strcmp(up.paramSet.ref_method, 'rrs')
        rr_ref = find_rr_ref_from_ref_rrs(wins, data, subj);
        save_or_append_data
        continue
    end
    
    %% Find ref RRs from simultaneous pressure respiratory signal
    if strcmp(up.paramSet.ref_method, 'paw')
        % Extract respiratory signals
        ext_sig = extract_resp_sigs(data, subj, up);
        % Determine SNRs
        snr_sig = calc_snrs(ext_sig, subj, wins, up);
        clear ext_sig
        % LPF to exclude freqs above resp
        lpf_sig = lpf_to_exclude_resp(snr_sig, subj, up);
        clear snr_sig
        % setup
        no_wins = length(wins.t_start);
        temp_t = mean([wins.t_start(:)' ; wins.t_end(:)']); temp_t = temp_t(:);
        rr_ref.v = nan(length(wins.t_start),1);
        rr_ref.snr_log = false(length(wins.t_start),1);
        % cycle through each window
        breath_times = [];
        for win_no = 1 : no_wins
                % select data for this window
                rel_els = find(lpf_sig.t >= wins.t_start(win_no) & lpf_sig.t <= wins.t_end(win_no));
                rel_data.t = lpf_sig.t(rel_els);
                rel_data.v = lpf_sig.v(rel_els);
                % interpolate data to be at fixed time values
                downsample_freq = 5;
                interp_data.t = wins.t_start(win_no): (1/downsample_freq) : wins.t_end(win_no);
                interp_data.v = interp1(rel_data.t, rel_data.v, interp_data.t, 'linear');
                % get rid of nans
                interp_data.v(isnan(interp_data.v)) = median(interp_data.v(~isnan(interp_data.v)));
                % normalise data
                rel_sig.v = (interp_data.v - mean(interp_data.v))/std(interp_data.v);
                rel_sig.t = interp_data.t;
                rel_sig.snr = lpf_sig.snr;
                rel_sig.fs = downsample_freq;
                clear rel_data interp_data
                % Identify positive gradient threshold-crossings
                [rr_ref.v(win_no), rr_ref.snr_log(win_no), temp_breath_times] = pos_grad_thresh(rel_sig, wins, win_no, up);
                rr_ref.t = temp_t;
                % store breath times
                breath_times =[breath_times; temp_breath_times(:)];
                
                clear rel_sig ave_breath_duration win_breaths el sig downsample_freq rel_els rel_paw rel_imp
                    
%                 if ~isempty(strfind(up.paths.data_load_filename, 'LISTEN'))
%                     
%                     %% CtO
%                     pic_save_name = [data(subj).subj, '_' data(subj).group, '_win_' num2str(win_no)];
%                     [qual, rr_cto] = ref_cto_mod(rel_sum_both, up, pic_save_name);
%                     
%                     %% Final RR for this window
%                     rr_ref.comb_imp.v(win_no,1) = rr_cto;
%                     rr_ref.comb_imp.snr_log(win_no,1) = qual;
%                     
%                     clear qual pic_save_name rr_cto rr_fft rrs breaths breath_times sum_both ave_breath_duration win_breaths el sig rel_data interp_data downsample_freq rel_els rel_paw rel_imp
%                     
%                 else
%                     
%                     %% CtO
%                     rr_cto = ref_cto(rel_sum_both, up);
%                     
%                     %% FFT
%                     rr_fft = ref_fft(rel_sum_both, downsample_freq, up);
%                     rr_wch = ref_wch(rel_sum_both, downsample_freq, up);
%                     
%                     %% Final RR for this window
%                     
%                     rrs = [rr_cto, rr_fft];
%                     if range(rrs) < 2
%                         rr_ref.comb_imp.v(win_no,1) = mean(rrs);
%                         rr_ref.comb_imp.snr_log(win_no,1) = true;
%                     else
%                         rr_ref.comb_imp.v(win_no,1) = nan;
%                         rr_ref.comb_imp.snr_log(win_no,1) = false;
%                     end
%                     
%                     clear rr_cto rr_fft rrs breaths breath_times sum_both ave_breath_duration win_breaths el sig rel_data interp_data downsample_freq rel_els rel_paw rel_imp
%                 end
                
        end
        clear no_wins s wins win_no temp_t lpf_sig
        
%         if ~isempty(strfind(up.paths.data_load_filename, 'VORTAL'))
%             for inputs = {'sum_both', 'imp', 'paw'}
%                 eval(['rel_data = rr_ref.' inputs{1,1} '_thresh;'])
%                 if sum(strcmp(fieldnames(rel_data), 'snr')) && strcmp(inputs{1,1}, 'imp')
%                     rel_data.resp_sqi = rel_data.snr;
%                     rel_data.snr_log = logical(rel_data.snr > 0.75);
%                 elseif sum(strcmp(fieldnames(rel_data), 'snr'))
%                     rel_data.snr_log = logical(rel_data.snr);
%                     rel_data = rmfield(rel_data, 'snr');
%                 end
%                 eval(['rr_ref.' inputs{1,1} '_thresh = rel_data;'])
%                 clear rel_data
%             end
%         end
        
        %% Save ref RRs to file
        save_or_append_data
        master_breath_times{subj} = breath_times;
        clear rr_ref breath_times
    end
    
    %% Find ref RRs from simultaneous impedance respiratory signal
    if strcmp(up.paramSet.ref_method, 'imp')
        % Extract respiratory signals
        ext_sig = extract_resp_sigs(data, subj, up);
        % Determine SNRs
        snr_sig = calc_snrs(ext_sig, subj, wins, up);
        clear ext_sig
        % LPF to exclude freqs above resp
        lpf_sig = lpf_to_exclude_resp(snr_sig, subj, up);
        clear snr_sig
        % setup
        no_wins = length(wins.t_start);
        temp_t = mean([wins.t_start(:)' ; wins.t_end(:)']); temp_t = temp_t(:);
        rr_ref.v = nan(length(wins.t_start),1);
        rr_ref.snr_log = false(length(wins.t_start),1);
        % cycle through each window
        for win_no = 1 : no_wins
            % select data for this window
            rel_els = find(lpf_sig.t >= wins.t_start(win_no) & lpf_sig.t <= wins.t_end(win_no));
            rel_data.t = lpf_sig.t(rel_els);
            rel_data.v = lpf_sig.v(rel_els);
            % interpolate data to be at fixed time values
            downsample_freq = 5;
            interp_data.t = wins.t_start(win_no): (1/downsample_freq) : wins.t_end(win_no);
            interp_data.v = interp1(rel_data.t, rel_data.v, interp_data.t, 'linear');
            % get rid of nans
            interp_data.v(isnan(interp_data.v)) = median(interp_data.v(~isnan(interp_data.v)));
            % normalise data
            rel_sig.v = (interp_data.v - mean(interp_data.v))/std(interp_data.v);
            rel_sig.t = interp_data.t;
            rel_sig.snr = lpf_sig.snr;
            rel_sig.fs = downsample_freq;
            clear rel_data interp_data
            % perform combined CtO and FFT analysis:
            % CtO
            rr_cto = ref_cto(rel_sig, up);
            % FFT
            rr_fft = ref_fft(rel_sig, downsample_freq, up);
            % Final RR for this window
            rrs = [rr_cto, rr_fft];
            if range(rrs) < 2
                rr_ref.v(win_no,1) = mean(rrs);
                rr_ref.snr_log(win_no,1) = true;
            else
                rr_ref.v(win_no,1) = nan;
                rr_ref.snr_log(win_no,1) = false;
            end
            
            clear rel_sig rr_cto rr_fft rrs rel_data interp_data downsample_freq rel_els
            
        end
        clear no_wins s wins win_no temp_t lpf_sig
        
%         if ~isempty(strfind(up.paths.data_load_filename, 'VORTAL'))
%             for inputs = {'sum_both', 'imp', 'paw'}
%                 eval(['rel_data = rr_ref.' inputs{1,1} '_thresh;'])
%                 if sum(strcmp(fieldnames(rel_data), 'snr')) && strcmp(inputs{1,1}, 'imp')
%                     rel_data.resp_sqi = rel_data.snr;
%                     rel_data.snr_log = logical(rel_data.snr > 0.75);
%                 elseif sum(strcmp(fieldnames(rel_data), 'snr'))
%                     rel_data.snr_log = logical(rel_data.snr);
%                     rel_data = rmfield(rel_data, 'snr');
%                 end
%                 eval(['rr_ref.' inputs{1,1} '_thresh = rel_data;'])
%                 clear rel_data
%             end
%         end
        
        %% Save ref RRs to file
        save_or_append_data
        clear rr_ref
    end
    
end

end

function rr_ref = find_rr_ref_from_ref_rrs(wins, data, subj)
rr_ref.t = mean([wins.t_start(:)' ; wins.t_end(:)']); rr_ref.t = rr_ref.t(:);
rr_ref.v = nan(length(rr_ref.t),1);
rel_data = data(subj).ref.params.rr;
for win_no = 1 : length(wins.t_start)
    rel_els = rel_data.t >= wins.t_start(win_no) & rel_data.t < wins.t_end(win_no);
    rr_ref.v(win_no) = nanmean(rel_data.v(rel_els));
end
end

function rr_ref = find_rr_ref_from_ref_breaths(wins, data, subj, up)
rel_breath_timings = data(subj).reference.breaths.t;

rr_ref.t = mean([wins.t_start(:)' ; wins.t_end(:)']); rr_ref.t = rr_ref.t(:);
rr_ref.v = nan(length(wins.t_start),1);
for win_no = 1 : length(rr_ref.t)
    win_breaths = rel_breath_timings >= wins.t_start(win_no) ...
        & rel_breath_timings < wins.t_end(win_no);
    if sum(win_breaths) == 0
        rr_ref.v(win_no) = NaN;
    else
        ave_breath_duration = range(rel_breath_timings(win_breaths))/(sum(win_breaths)-1);
        rr_ref.v(win_no) = 60/ave_breath_duration;
        % rr.v(win_no) = mean(60./diff(rr_breaths));   % doesn't work as outliers have a big effect (eg spurious breath detections)
    end
end

% in the absence of snrs, assume they are all good:
ideal_snr = 10;
rr_ref_snr = ideal_snr*ones(length(rr_ref.v),1);
rr_ref.snr_log = logical(rr_ref_snr > up.paramSet.paw_snr_thresh);

% mark those which have implausible RRs as bad:
bad_els = rr_ref.v>up.paramSet.rr_range(2) | rr_ref.v<up.paramSet.rr_range(1);
rr_ref.snr_log(bad_els) = false;
end

function new_sig = extract_resp_sigs(data, subj, up)

%% See which resp signals are present:
if sum(strcmp(up.paramSet.ref_method, 'imp'))
    % select impedance:
    imp.v = data(subj).ref.resp_sig.imp.v;
    imp.fs = data(subj).ref.resp_sig.imp.fs;
    imp.t = (1/imp.fs)*(1:length(imp.v));
    new_sig = imp;
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
lpf_sig.v = lp_filter_signal_to_remove_freqs_above_resp(d_s_win.v, d_s_win.fs, up);
lpf_sig.fs = d_s_win.fs;
lpf_sig.snr = sig.snr;
% NB: if you use the e_vlf signal then the freqs below resp have already been removed.

end

function [rr_ref_val, rr_ref_snr, breath_times] = pos_grad_thresh(rel_sig, wins, win_no, up, sig_name)

% Identify peak detection threshold
if strcmp(up.paramSet.ref_method, 'imp')
    thresh = up.paramSet.resp_sig_thresh.imp;
    rr_ref_snr = true;
elseif strcmp(up.paramSet.ref_method, 'paw')
    thresh = up.paramSet.resp_sig_thresh.paw;
    rr_ref_snr = logical(rel_sig.snr(win_no) > up.paramSet.paw_snr_thresh);
end

% Find breath times
breath_times = [];
for el = 1 : (length(rel_sig.v)-1)
    if rel_sig.v(el) < thresh & rel_sig.v(el+1) > thresh
        breath_times = [breath_times, mean([rel_sig.t(el), rel_sig.t(el+1)])];
    end
end

% Find RR
win_breaths = breath_times >= wins.t_start(win_no) ...
    & breath_times < wins.t_end(win_no);
if sum(win_breaths) == 0
    rr_ref_val = NaN;
else
    ave_breath_duration = range(breath_times(win_breaths))/(sum(win_breaths)-1);
    rr_ref_val = 60/ave_breath_duration;
end

end

function rr_cto = ref_cto(sum_both, up)

% identify peaks
diffs_on_left_of_pt = diff(sum_both.v); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt>0);
diffs_on_right_of_pt = diff(sum_both.v); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt<0);
peaks = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
% identify troughs
diffs_on_left_of_pt = diff(sum_both.v); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt<0);
diffs_on_right_of_pt = diff(sum_both.v); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt>0);
troughs = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
% define threshold
q3 = quantile(sum_both.v(peaks), 0.75);
thresh = 0.2*q3;
% find relevant peaks and troughs
extrema = sort([peaks(:); troughs(:)]);
rel_peaks = peaks(sum_both.v(peaks) > thresh);
rel_troughs = troughs(sum_both.v(troughs) < 0);

% find valid breathing cycles
% valid cycles start with a peak:
valid_cycles = zeros(length(rel_peaks)-1,1);
cycle_durations = nan(length(rel_peaks)-1,1);
for peak_no = 1 : (length(rel_peaks)-1)
    
    % valid if there is only one rel_trough between this peak and the
    % next
    cycle_rel_troughs = rel_troughs(rel_troughs > rel_peaks(peak_no) & rel_troughs < rel_peaks(peak_no+1));
    if length(cycle_rel_troughs) == 1
        valid_cycles(peak_no) = 1;
        cycle_durations(peak_no) = sum_both.t(rel_peaks(peak_no+1)) - sum_both.t(rel_peaks(peak_no));
    end
    
end

% Calc RR
if sum(valid_cycles) == 0
    rr_cto = nan;
else
    % Using average breath length
    ave_breath_duration = nanmean(cycle_durations);
    rr_cto = 60/ave_breath_duration;
end

end

function rr_zex = ref_zex(sum_both, up)

%% identify individual breaths from the raw signal using +ve grad zero-crossing detection:
val_on_left_of_pt = sum_both.v(1:(end-1)); left_log = logical(val_on_left_of_pt < 0);
val_of_pt = sum_both.v(2:end); pt_log = logical(val_of_pt >= 0);
breaths = find(left_log & pt_log)+1;
if isempty(breaths)
    rr_zex = nan;
    return
end

%% Calc RR
% Only using time period spanning a whole number of breaths
win_length = sum_both.t(breaths(end)) - sum_both.t(breaths(1));
rr_zex = 60*(length(breaths)-1)/win_length;

end

function rr_wch = ref_wch(sum_both, downsample_freq, up)

segLen = 2^nextpow2(12*downsample_freq);
noverlap = segLen/2;
[data.power, data.freqs] = pwelch(sum_both.v,segLen,noverlap, [], downsample_freq);

% Find spectral peak
[rr_wch, ~, ~] = find_spectral_peak(data, up);

end

function rr_fft = ref_fft(sum_both, downsample_freq, up)

% Find FFT
WINLENGTH = length(sum_both.v);
NFFT = 2^nextpow2(WINLENGTH);
HAMMWIN = hamming(WINLENGTH);
HAMMWIN = HAMMWIN(:);
f_nyq = downsample_freq/2;
FREQS = f_nyq.*linspace(0, 1, NFFT/2+1);            % Array of correspondent FFT bin frequencies, in BR (RPM)
WINDATA = detrend(sum_both.v(:));                      % Remove the LSE straight line from the data
WINDATA = WINDATA .* HAMMWIN;
myFFT = fft(WINDATA, NFFT);
myFFT = myFFT(1 : NFFT/2 + 1);
myFFT = 2.*abs(myFFT/NFFT);
psdx = (1/(downsample_freq*NFFT)) * abs(myFFT).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
power = 10*log10(psdx); power = power(:);
freqs = FREQS; freqs = freqs(:);

% Find respiratory peak
freq_range = up.paramSet.rr_range/60;
cand_els = zeros(length(power),1);
for s = 2 : (length(power)-1)
    if power(s) > power(s-1) & power(s) > power(s+1) & freqs(s) > freq_range(1) & freqs(s) < freq_range(2)
        cand_els(s) = 1;
    end
end
cand_els = find(cand_els);

[~, r_el] = max(power(cand_els));
r_el = cand_els(r_el);
r_freq = freqs(r_el);
if ~isempty(r_freq)
    rr_fft = 60*r_freq;
else
    rr_fft = nan;
end

end

function calc_spec(sig, downsample_freq, up)

close all

% Find FFT
WINLENGTH = length(sig.v);
NFFT = 2^nextpow2(WINLENGTH);
HAMMWIN = hamming(WINLENGTH);
HAMMWIN = HAMMWIN(:);
f_nyq = downsample_freq/2;
FREQS = f_nyq.*linspace(0, 1, NFFT/2+1);            % Array of correspondent FFT bin frequencies, in BR (RPM)
WINDATA = detrend(sig.v(:));                      % Remove the LSE straight line from the data
WINDATA = WINDATA .* HAMMWIN;
myFFT = fft(WINDATA, NFFT);
myFFT = myFFT(1 : NFFT/2 + 1);
myFFT = 2.*abs(myFFT/NFFT);
psdx = (1/(downsample_freq*NFFT)) * abs(myFFT).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
power = 10*log10(psdx); power = power(:);
freqs = FREQS; freqs = freqs(:);

plot(freqs, power)

end

function PrintFigs(h, paper_size, savepath, up)
set(h,'PaperUnits','inches');
set(h,'PaperSize', [paper_size(1), paper_size(2)]);
set(h,'PaperPosition',[0 0 paper_size(1) paper_size(2)]);
% print(h,'-dpdf',savepath)
print(h,'-dpng',savepath)

% % you need to download 'export_fig' from:
% % http://uk.mathworks.com/matlabcentral/fileexchange/23629-export-fig
% export_fig_dir_path = 'C:\Users\pc13\Documents\GitHub\phd\Tools\Other Scripts\export_fig\';
% addpath(export_fig_dir_path)
% export_fig(savepath, '-eps')

close all
end