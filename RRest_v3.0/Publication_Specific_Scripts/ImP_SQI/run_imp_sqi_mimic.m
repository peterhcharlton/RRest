function run_imp_sqi_mimic
% run_imp_sqi_mimic analyses the performance of the impedance pneumography
% signal quality index (SQI) using data from the MIMIC database.
%
%               run_imp_sqi_mimic
%
%	Inputs:
%       mimic_imp_sqi_data.mat - a file containing curated and
%          annotated data from the MIMIC database.
%       Specify the path to this file in the "setup_universal_params"
%       function below.
%
%	Outputs:
%       The script outputs the numerical results reported in the
%       publication verbatim (in the Command Window). In addition, figures
%       are created (including those provided in the publication).
%           
%   Further Information:
%       This version of run_imp_sqi_mimic is provided to facilitate
%       reproduction of the analyses reported in:
%           Charlton P. H. et al., "An impedance pneumography signal
%           quality index for respiratory rate monitoring: design,
%           assessment and application", [under review]
%       Further information on this study can be obtained at:
%           http://peterhcharlton.github.io/RRest/imp_sqi.html
%
%   Comments, Questions, Criticisms, Feedback, Contributions:
%       See: http://peterhcharlton.github.io/RRest/contributions.html
%
%   Version:
%       v.0.1 - accompanying peer review, 6th Aug 2020 by Peter H Charlton
%
%   Licence:
%       This program is available under the GNU public license. This
%       program is free software: you can redistribute it and/or modify 
%       it under the terms of the GNU General Public License as published by
%       the Free Software Foundation, either version 3 of the License, or
%       (at your option) any later version. This program is distributed in
%       the hope that it will be useful, but WITHOUT ANY WARRANTY; without
%       even the implied warranty of MERCHANTABILITY or FITNESS FOR A
%       PARTICULAR PURPOSE.  See the GNU General Public License for more
%       details: <http://www.gnu.org/licenses/>.
%

%% Setup
% Declare Universal Parameters (includes the path of the collated data file)
up = setup_universal_params;

%% Load Data 
% Load data containing impedance signal and reference breath annotations
load_data(up);

%% SQI Evaluation
% Determine the performance of each algorithm
AnalysePerformances(up)

end

function up = setup_universal_params

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% USER-SPECIFIC SETTINGS

% Specify path of folder containing the collated and annotated dataset
%   ("mimic_imp_sqi_data.mat")
up.paths.root_folder = '/Users/petercharlton/Documents/Data/ImP_SQI_mimic_dataset_test/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compulsory setup tasks

% set current directory to that of this file
cd(fileparts(mfilename('fullpath')))
% Add all functions corresponding to options for each component to the path
addpath(genpath(fileparts(mfilename('fullpath'))));
% Specify name of universal params file:
up_filename = 'up';
up_path = [fileparts(mfilename('fullpath')), filesep, up_filename];
% Specify whether or not to re-make the universal_params
redo_up = 1;

%% Optional setup tasks

% Don't remake the universal_params file unless needed:
if exist([up_filename, '.mat'], 'file') & ~redo_up
    fprintf('\n--- Loading Universal Parameters ');
    load(up_filename);
    return
end

% Otherwise make the universal_params file:
fprintf('\n--- Creating Universal Parameters ');

%% Specify the folders in which to save results files

% File paths
up.paths.raw_imp_data = [up.paths.root_folder, 'mimic_imp_sqi_data.mat'];
up.paths.data_save_folder = up.paths.root_folder;
up.paths.plots_save_folder = [up.paths.data_save_folder, 'Figures', filesep];
up.paths.filenames.data_mat = 'data_mat';
up.paths.filenames.win_timings = '_wins';
up.paths.filenames.analysis = 'analysis';
up.paths.filenames.algorithms = 'algorithms';
up.paths.filenames.evaluation = 'eval_data';

% Make directories
directories_to_make = {up.paths.data_save_folder, up.paths.plots_save_folder};
for s = 1 : length(directories_to_make)
    rel_path = directories_to_make{s};
    if ~exist(rel_path, 'dir')
        mkdir(rel_path)
    end    
end

% create subject list
load(up.paths.raw_imp_data);
up.paramSet.subj_list = 1 : length(data); clear data

% window parameters
up.paramSet.winLeng = 32;     % the duration of each window in secs
up.paramSet.winStep = 0; %up.paramSet.winLeng/2; %3;      % the number of secs between each consecutive window

% Eliminate HFs (above resp freqs)
up.paramSet.elim_hf.Fpass = 1.2;  % in Hz
up.paramSet.elim_hf.Fstop = 0.895;  % in Hz
up.paramSet.elim_hf.Dpass = 0.057501127785;
up.paramSet.elim_hf.Dstop = 0.01;

% duration of Tukey window taper in secs
up.paramSet.tukey_win_duration_taper = 2;

% Range of plausible RR freqs
up.paramSet.rr_range = [4 60];

% Ref RR threshold
up.paramSet.resp_sig_thresh.paw = 0.42; % thresh for use of paw (just a single paw, not 2*paw)
up.paramSet.resp_sig_thresh.imp = 0.1; % thresh for use of imp  (just a single imp, not 2*imp)
% Paw SNR
up.paramSet.paw_snr_thresh = -8;

% Evaluation
up.eval.rel_thresh.Gr = linspace(0.95, 1.00, 100);
up.eval.rel_thresh.Qw = linspace(0.90, 1.00, 100);
up.eval.rel_thresh.QIr = linspace(0.5, 3.00, 100);
up.eval.skip_plots = 0;

%% save universal params file
save(up_path, 'up');

end

function load_data(up)
fprintf('\n--- Loading Data, Calculating RRs and Calculating Novel and Agreement SQIs');

%% Load data
load(up.paths.raw_imp_data);

%% Cycle through each subject
counter_no = 0;
[res.subj, res.flat_line, res.rr_ann, res.sqi_novel, res.rr_novel, res.rr_agree, res.rr_mon, res.rr_wch, res.rr_cto, res.prop_norm_dur, res.prop_bad_breaths, res.R2, res.R2min, res.win_no] = deal(nan(20000,1));
[res.id] = deal(cell(20000,1));
res.for_rr_perf_analysis = false(20000,1);
for s = up.paramSet.subj_list
    
    %% See if it has been annotated
    if sum(strcmp(fieldnames(data(s).ref), 'ann'))
        curr_ann = data(s).ref.ann;
    else
        curr_ann = struct;
    end
    
    %% Extract impedance signal for this subject
    imp.fs = data(s).ref.resp_sig.imp.fs;
    imp.v = data(s).ref.resp_sig.imp.v;
    imp.t = [0:(length(imp.v)-1)]/imp.fs;
    monitor_rr.v = double(data(s).ref.params.rr.v);
    monitor_rr.t = double(data(s).ref.params.rr.t);
    
    %% Filter impedance signal
    % Filtering impedance signal to eliminate high frequencies above respiration
    imp_filt = lpf_to_exclude_resp(imp, up);
    
    %% Make plot of subject's data
    do_plot = 0;
    if do_plot
        
        figure('Position', [20,20,800,500])
        ftsize = 20; lwidth = 2;
        subplot(2,1,1), plot(imp.t, imp.v),
        ylab = ylabel('ImP', 'FontSize', ftsize, 'Rotation', 0);
        set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
        set(gca, 'FontSize', ftsize)
        xlim([0 10])
        subplot(2,1,2),
        plot(monitor_rr.t, monitor_rr.v, 'LineWidth', lwidth),
        ylab = ylabel({'RR', 'bpm'}, 'FontSize', ftsize, 'Rotation', 0);
        set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
        xlabel('Time [s]', 'FontSize', ftsize)
        set(gca, 'FontSize', ftsize)
        xlim([0 10])
        ylim([min(monitor_rr.v)-0.1*range(monitor_rr.v), max(monitor_rr.v)+0.1*range(monitor_rr.v)])
        
        savepath = [up.paths.plots_save_folder, 'sample_recording'];
        PrintFigs(gcf, paper_size, savepath, up)
        close all
    end
    clear do_plot
    
    %% Normalise signals in each window
    % The impedance signal corresponding to each window is extracted, then
    % it is normalised to have a mean of 0 and standard deviation of 1.
    win_starts = imp.t(1):up.paramSet.winLeng:imp.t(end);
    for win_no = 1 : length(win_starts)
        
        counter_no = counter_no + 1;
        
        % input subject
        res.win_no(counter_no) = win_no;
        res.subj(counter_no) = s;
        res.id{counter_no} = data(s).fix.id;
        
        % select data for this window
        win_start = win_starts(win_no);
        win_end = win_start + up.paramSet.winLeng;
        rel_els = find(imp_filt.t >= win_start & imp_filt.t <= win_end);
        rel_data.t = imp_filt.t(rel_els);
        rel_data.v = imp_filt.v(rel_els);
        
        % identify any annotations within this window:
        if sum(strcmp(fieldnames(curr_ann), 'breath_els'))
            rel_ann_els = curr_ann.breath_els(curr_ann.breath_els >= rel_els(1) & curr_ann.breath_els <= rel_els(end));
            if sum(rel_ann_els) > 1
                first_breath = min(imp.t(rel_ann_els));
                last_breath = max(imp.t(rel_ann_els));
                no_breaths = length(rel_ann_els)-1;
                res.rr_ann(counter_no) = 60/((last_breath-first_breath)/no_breaths);
            else
                res.rr_ann(counter_no) = nan;
            end
            clear rel_ann_els
        else
            res.rr_ann(counter_no) = nan;
        end
        
        % say whether or not this is in the first five annotated for this subject 
        subj_els = res.subj == s;
        if sum(~isnan(res.rr_ann(subj_els)))<=5 && ~isnan(res.rr_ann(counter_no))
            res.for_rr_perf_analysis(counter_no) = true;
        end
        
        % determine whether it was a flat line - if so then skip this window
        if length(unique(imp.v(rel_els))) == 1 || sum(isnan(imp.v(rel_els)))>=1
            res.flat_line(counter_no) = 1;
            continue
        end
        res.flat_line(counter_no) = 0;
        % interpolate data to be at fixed time values
        downsample_freq = 5;
        interp_data.t = win_start: (1/downsample_freq) : win_end;
        interp_data.v = interp1(rel_data.t, rel_data.v, interp_data.t, 'linear');
        % get rid of nans
        interp_data.v(isnan(interp_data.v)) = median(interp_data.v(~isnan(interp_data.v)));
        % normalise data
        rel_imp.v_n = (interp_data.v - mean(interp_data.v))/std(interp_data.v);
        rel_imp.t_n = interp_data.t;
        rel_imp.fs = downsample_freq;
        clear interp_data rel_data downsample_freq rel_els
        
        %% Estimate RR from the impedance signal
        
        % - using the novel modified count-orig method (which also provides signal quality metrics)
        [res.sqi_novel(counter_no), res.rr_novel(counter_no), res.prop_norm_dur(counter_no), res.prop_bad_breaths(counter_no), res.R2(counter_no), res.R2min(counter_no)] ...
            = ref_cto_mod(rel_imp, up, 'no');
        
        % - using the original count-orig method (returns a nan if poor-quality)
        res.rr_cto(counter_no,1) = ref_cto(rel_imp, up);
        
        % - using an FFT
        res.rr_wch(counter_no,1) = ref_wch(rel_imp, up);
        
        % using the agreement SQI
        res.rr_agree(counter_no,1) = mean([res.rr_cto(counter_no), res.rr_wch(counter_no)]);
        
        %% Estimate RR from monitor numerics
        
        rel_els = find(monitor_rr.t >= win_start & monitor_rr.t <= win_end);
        rel_data.t = monitor_rr.t(rel_els);
        rel_data.v = monitor_rr.v(rel_els);
%         res.rr_mon(counter_no,1) = nanmean(rel_data.v);
        res.rr_mon(counter_no,1) = nanmedian(rel_data.v);
        
        clear rel_els rel_data win_start win_end  rel_imp
        
    end    
    clear imp imp_filt no_wins qual monitor_rr curr_ann
    
end

% Eliminate nans
fields = fieldnames(res);
rel_rows = 1:counter_no-1;
for field_no = 1 : length(fields)
    eval(['res.' fields{field_no} ' = res.' fields{field_no} '(rel_rows);'])    
end
clear s counter_no field_no fields rel_rows temp


%% Calculate agreement SQI
fprintf('\n--- Calculating agreement SQI ');
good_els = abs(res.rr_cto-res.rr_wch) < 2;
res.agree_sqi = false(length(res.rr_cto),1);
res.agree_sqi(good_els) = true;
clear good_els bad_els

for period = {'all'} %{'training', 'validation', 'all'}
    if strcmp(period{1,1}, 'training')
        rel_rows = res.subj <= ceil(max(res.subj)/2);
    elseif strcmp(period{1,1}, 'validation')
        rel_rows = res.subj > ceil(max(res.subj)/2);
    else
        rel_rows = ~isnan(res.subj);
    end
    temp.design_header = {'R2', 'R2min', 'prop_norm_dur', 'prop_bad_breaths'};
    temp.design_data = [res.R2(rel_rows,:), res.R2min(rel_rows,:), res.prop_norm_dur(rel_rows,:), res.prop_bad_breaths(rel_rows,:)];
    temp.response_header = {'subj', 'flat_line', 'rr_ann', 'rr_novel', 'sqi_novel', 'sqi_agree', 'rr_agree', 'rr_mon', 'win_no', 'for_rr_perf_analysis'};
    temp.response_data = [res.subj(rel_rows,:), res.flat_line(rel_rows,:), res.rr_ann(rel_rows,:), res.rr_novel(rel_rows,:), res.sqi_novel(rel_rows,:), res.agree_sqi(rel_rows,:), res.rr_agree(rel_rows,:), res.rr_mon(rel_rows,:), res.win_no(rel_rows,:), res.for_rr_perf_analysis(rel_rows,:)];
    eval(['data_mat.' period{1,1} ' = temp;']);
    clear temp
end

% Find out percentage of segments which are flat lines:
rel_col = strcmp(data_mat.all.response_header, 'flat_line');
perc_flat_line = 100*sum(data_mat.all.response_data(:,rel_col) == 1)/length(data_mat.all.response_data);
fprintf(['\n - Percentage of flat line segments: ' num2str(perc_flat_line) ' %% \n']);

% Find out number of annotated segments
rel_col = strcmp(data_mat.all.response_header, 'rr_ann');
ann_segs = sum(~isnan(data_mat.all.response_data(:,rel_col)));
fprintf(['\n - Number of annotated segments: ' num2str(ann_segs) '\n']);

%% Save data to file
save_name = up.paths.filenames.data_mat;
savepath = [up.paths.data_save_folder, up.paths.filenames.data_mat, '.mat'];
save(savepath, 'data_mat');

end

function imp = lpf_to_exclude_resp(imp, up)

imp.v(isnan(imp.v)) = median(imp.v(~isnan(imp.v)));

%% Window signal to reduce edge effects
duration_of_signal = imp.t(end) - imp.t(1);
prop_of_win_in_outer_regions = 2*up.paramSet.tukey_win_duration_taper/duration_of_signal;
tukey_win = tukeywin(length(imp.v), prop_of_win_in_outer_regions);
d_s_win = imp;    % copy time and fs
d_s_win.v = detrend(imp.v(:)).*tukey_win(:);

%% LPF to remove freqs above resp
respWave.t = d_s_win.t;
respWave.v = lp_filter_signal_to_remove_freqs_above_resp(d_s_win.v, d_s_win.fs, up);
respWave.fs = d_s_win.fs;
imp = respWave;
end

function s_filt = lp_filter_signal_to_remove_freqs_above_resp(s, Fs, up)
%% Filter pre-processed signal to remove freqs above resp

% parameters for the low-pass filter to be used
flag  = 'scale';
Dpass = up.paramSet.elim_hf.Dpass;
Dstop = up.paramSet.elim_hf.Dstop;
Fstop = up.paramSet.elim_hf.Fstop;
Fpass = up.paramSet.elim_hf.Fpass;

% create filter
[N,Wn,BETA,TYPE] = kaiserord([Fstop Fpass]/(Fs/2), [1 0], [Dstop Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.0159;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(Fs/2);

% Prepare signal
s_dt=detrend(s);
s_filt = filtfilt(AMfilter.numerator, 1, s_dt);
end

function [qual, rr_cto, prop_norm_dur, prop_bad_breaths, R2, R2min] = ref_cto_mod(sum_both, up, save_name)

sum_both.t_n = sum_both.t_n-sum_both.t_n(1);
sum_both.v_n = -1*detrend(sum_both.v_n);

%% Identify relevant peaks and troughs

% identify peaks
diffs_on_left_of_pt = diff(sum_both.v_n); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt>0);
diffs_on_right_of_pt = diff(sum_both.v_n); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt<0);
peaks = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
% identify troughs
diffs_on_left_of_pt = diff(sum_both.v_n); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt<0);
diffs_on_right_of_pt = diff(sum_both.v_n); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt>0);
troughs = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
% define peaks threshold
q3 = quantile(sum_both.v_n(peaks), 0.75);
thresh = 0.2*q3;
% find relevant peaks
rel_peaks = peaks(sum_both.v_n(peaks) > thresh);
% define troughs threshold
q3t = quantile(sum_both.v_n(troughs), 0.25);
thresh = 0.2*q3t;
% find relevant troughs
rel_troughs = troughs(sum_both.v_n(troughs) < thresh);

%% find valid breathing cycles
% exclude peaks which aren't the highest between a pair of consecutive
% troughs:
invalid_peaks = zeros(length(rel_peaks),1);
for trough_pair_no = 1 : (length(rel_troughs)-1)
    
    % identify peaks between this pair of troughs
    cycle_rel_peak_els = find(rel_peaks > rel_troughs(trough_pair_no) & rel_peaks < rel_troughs(trough_pair_no+1));
    cycle_rel_peaks = rel_peaks(cycle_rel_peak_els);
    if length(cycle_rel_peaks) > 1
        [~, rel_el] = max(sum_both.v_n(cycle_rel_peaks));
        bad_rel_peaks_els = setxor(1:length(cycle_rel_peak_els), rel_el);
        invalid_peaks(cycle_rel_peak_els(bad_rel_peaks_els)) = 1;
    end
end
rel_peaks = rel_peaks(~invalid_peaks);

% if there is more than one initial peak (i.e. before the first trough) then take the highest:
initial_peaks = find(rel_peaks < rel_troughs(1));
other_peaks = find(rel_peaks >= rel_troughs(1));
if length(initial_peaks)>1
    [~, rel_initial_peak] = max(sum_both.v_n(rel_peaks(initial_peaks)));
    rel_peaks = rel_peaks([rel_initial_peak, other_peaks]);
end

% valid cycles start with a peak:
valid_cycles = false(length(rel_peaks)-1,1);
cycle_durations = nan(length(rel_peaks)-1,1);

for peak_no = 2 : length(rel_peaks)
    
    % exclude if there isn't a rel trough between this peak and the
    % previous one
    cycle_rel_troughs = rel_troughs(rel_troughs > rel_peaks(peak_no-1) & rel_troughs < rel_peaks(peak_no));
    if length(cycle_rel_troughs) ~= 0
        valid_cycles(peak_no-1) = true;
        cycle_durations(peak_no-1) = sum_both.t_n(rel_peaks(peak_no)) - sum_both.t_n(rel_peaks(peak_no-1));
    end
end
valid_cycle_durations = cycle_durations(valid_cycles);

% Calc RR
if isempty(valid_cycle_durations)
    rr_cto = nan;
else
    % Using average breath length
    ave_breath_duration = mean(valid_cycle_durations);
    rr_cto = 60/ave_breath_duration;
end

%% Resiratory SQI

if isnan(rr_cto)
    qual = false;
    prop_norm_dur = 0;
    prop_bad_breaths = 100;
    R2 = 0;
    R2min = 0;
else
    
    %find mean breath-to-breath interval to define size of template
    rr=floor(mean(diff(rel_peaks)));
    ts=[];
    j=find(rel_peaks>rr/2);
    l=find(rel_peaks+floor(rr/2)<length(sum_both.v_n));
    new_rel_peaks = rel_peaks(j(1):l(end));
    if isempty(new_rel_peaks)
        qual = false;
        prop_norm_dur = 0;
        prop_bad_breaths = 100;
        R2 = 0;
        R2min = 0;
        return
    else
        %find breaths
        for k=1:length(new_rel_peaks)
            t=sum_both.v_n(new_rel_peaks(k)-floor(rr/2):new_rel_peaks(k)+floor(rr/2));
            tt=t/norm(t); tt = tt(:)';
            ts=[ts;tt];
        end
    end
    
    %find ave template
    if size(ts,1) > 1
        avtempl=mean(ts,1);
    else
        avtempl=nan(size(ts));
    end
    
    %now calculate correlation for every beat in this window
    r2 = nan(size(ts,1),1);
    for k=1:size(ts,1)
        r2(k)=corr2(avtempl,ts(k,:));
    end
    %calculate mean correlation coefficient
    R2=mean(r2);
    R2min = std(valid_cycle_durations)/mean(valid_cycle_durations);
    %peak_heights = sum_both.v_n(rel_troughs);
    %R2min = std(peak_heights)/mean(peak_heights);
    
    % calculate number of abnormal breath durations
    median_dur = median(valid_cycle_durations);
    temp = valid_cycle_durations > (1.5*median_dur) | valid_cycle_durations < (0.5*median_dur);
    prop_bad_breaths = 100*sum(temp)/length(temp);
    
    % find prop of window taken up by normal breath durations
    norm_dur = sum(valid_cycle_durations(~temp));
    win_length = sum_both.t_n(end) - sum_both.t_n(1);
    prop_norm_dur = 100*norm_dur/win_length;
    
    % determine whether this window is high or low quality
    if prop_norm_dur > 60 && prop_bad_breaths < 15 && R2 >= 0.75 && R2min < 0.25
        qual = true;
    else
        qual = false;
    end
    
end

%% Plot template and inidividual beats
save_name = 'sqi_process';
save_name = 'no';
if ~strcmp(save_name, 'no') && ~isnan(rr_cto)
    
    paper_size = [12, 8];
    figure('Position', [50, 50, 100*paper_size(1), 100*paper_size(2)], 'Color',[1 1 1])
    lwidth1 = 3; lwidth2 = 2; ftsize = 18;
    % plot signal
    subplot(2,2,[1,2]), plot(sum_both.t_n-sum_both.t_n(1), sum_both.v_n, 'LineWidth', lwidth2), hold on
    plot(sum_both.t_n(rel_peaks(logical([valid_cycles; 1])))-sum_both.t_n(1), sum_both.v_n(rel_peaks(logical([valid_cycles; 1]))), '.r', 'MarkerSize', 20)
    %plot(sum_both.t_n(new_rel_peaks)-sum_both.t_n(1), sum_both.v_n(new_rel_peaks), '.r', 'MarkerSize', 20)
    %plot(sum_both.t_n(rel_troughs)-sum_both.t_n(1), sum_both.v_n(rel_troughs), '.k', 'MarkerSize', 20)
    xlim([0, sum_both.t_n(end)-sum_both.t_n(1)])
    xlabel('Time [s]', 'FontSize', ftsize)
    ylab=ylabel('Imp', 'FontSize', ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    set(gca, 'FontSize', ftsize, 'YTick', [])
    % plot template
    time = 0:(length(avtempl)-1); time = time./sum_both.fs;
    subplot(2,2,3), hold on,
    for beat_no = 1 : size(ts,1)
        plot(time, ts(beat_no,:), 'color', 0.7*[1 1 1], 'LineWidth', lwidth2)
    end
    plot(time, avtempl, 'r', 'LineWidth', lwidth1)
    set(gca, 'YTick', [])
    xlabel('Time [s]', 'FontSize', ftsize)
    xlim([0, time(end)])
    ylab=ylabel('ImP', 'FontSize', ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    set(gca, 'FontSize', ftsize)
    %set ylim
    rang = range(ts(:));
    ylim([min(ts(:))-0.1*rang, max(ts(:))+0.1*rang]);
    % prop_norm_dur > 60 && prop_bad_breaths < 15 && R2 >= 0.75 && R2min < 0.25
    if R2 >= 0.75, col = 'b'; else, col = 'r'; end
    annotation('textbox',[0.5, 0.25, 0.1,0.1],'String',['R^2 = ' num2str(R2, 2)], 'Color', col, 'FontSize', ftsize, 'LineStyle', 'None')
    if prop_bad_breaths < 15, col = 'b'; else, col = 'r'; end
    annotation('textbox',[0.5, 0.18, 0.1,0.1],'String',['prop invalid breaths = ' num2str(prop_bad_breaths, 2), ' %'], 'Color', col, 'FontSize', ftsize, 'LineStyle', 'None')
    if prop_norm_dur > 60, col = 'b'; else, col = 'r'; end
    annotation('textbox',[0.5, 0.11, 0.1,0.1],'String',['prop valid duration = ' num2str(prop_norm_dur, 2), ' %'], 'Color', col, 'FontSize', ftsize, 'LineStyle', 'None')
    if R2min < 0.25, col = 'b'; else, col = 'r'; end
    annotation('textbox',[0.5, 0.04, 0.1,0.1],'String',['Breath interval variability = ' num2str(100*R2min, 2), ' %'], 'Color', col, 'FontSize', ftsize, 'LineStyle', 'None')
    
    %annotation('textbox',[0.5, 0.1, 0.1,0.1],'String',{['R^2 = ' num2str(R2, 2)] , ['prop breaths bad = ' num2str(prop_bad_breaths, 2) '%'], ['prop dur good = ' num2str(prop_norm_dur,2) '%'], ['norm SD durations = ' num2str(R2min,2) '%']}, 'FontSize', ftsize, 'LineStyle', 'None')
    if qual
        annotation('textbox',[0.8, 0.15, 0.1,0.1],'String','High Quality', 'Color', 'b', 'FontSize', ftsize+4, 'LineStyle', 'None')
    else
        annotation('textbox',[0.8, 0.15, 0.1,0.1],'String','Low Quality', 'Color', 'r', 'FontSize', ftsize+4, 'LineStyle', 'None')
    end
    savepath = [up.paths.plots_save_folder, save_name];
    PrintFigs(gcf, paper_size, savepath, up)
    close all
end

%% Plot summary figure
save_name = 'resp_sqi_summary_viva';
save_name = 'no';
if ~strcmp(save_name, 'no') && ~isnan(rr_cto)
    
    % check to see if a figure has been opened
    paper_size = [16, 7];
    if isempty(findall(0,'Type','Figure'))
        figure('Position', [50, 50, 100*paper_size(1), 100*paper_size(2)], 'Color',[1 1 1])
    end
    lwidth1 = 3; lwidth2 = 2; ftsize = 18; plotted = 0;
    
    % see if this signal is appropriate
    if R2 < 0.55
        % Plot signal
        subplot(2,3,1:2),plot(sum_both.t_n-sum_both.t_n(1), -1*sum_both.v_n, 'LineWidth', lwidth2), hold on
        plot(sum_both.t_n(rel_peaks(logical([valid_cycles; 1])))-sum_both.t_n(1), -1*sum_both.v_n(rel_peaks(logical([valid_cycles; 1]))), '.k', 'MarkerSize', 20), hold off
        plotted = 1;
        xlim([0, sum_both.t_n(end)-sum_both.t_n(1)]), ylab=ylabel('Imp', 'FontSize', ftsize, 'Rotation', 0);
        set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); set(gca, 'FontSize', ftsize, 'YTick', [])
        title('Low quality signal segment', 'FontSize', ftsize)
        %set ylim
        rang = range(sum_both.v_n);
        ylim(-1*fliplr([min(sum_both.v_n(:))-0.1*rang, max(sum_both.v_n(:))+0.1*rang]));
        % plot template
        subplot(2,3,3), cla, hold on
        time = 0:(length(avtempl)-1); time = time./sum_both.fs;
        for beat_no = 1 : size(ts,1)
            plot(time, -1*ts(beat_no,:), 'color', 0.7*[1 1 1], 'LineWidth', lwidth2)
        end
        plot(time, -1*avtempl, 'r', 'LineWidth', lwidth1), hold off
        set(gca, 'YTick', [])
        xlim([0, 6]) % Jan Submission
        xlim([0, time(end)])
        title(['Template (R^2 = ' num2str(R2,2) ')'], 'FontSize', ftsize)
        set(gca, 'FontSize', ftsize)
        %set ylim
        rang = range(ts(:));
        ylim(-1*fliplr([min(ts(:))-0.1*rang, max(ts(:))+0.1*rang]));
    elseif R2 >= 0.97
        % Plot signal
        subplot(2,3,4:5),plot(sum_both.t_n-sum_both.t_n(1), -1*sum_both.v_n, 'LineWidth', lwidth2), hold on
        plot(sum_both.t_n(rel_peaks(logical([valid_cycles; 1])))-sum_both.t_n(1), -1*sum_both.v_n(rel_peaks(logical([valid_cycles; 1]))), '.k', 'MarkerSize', 20), hold off
        plotted = 1;
        xlim([0, sum_both.t_n(end)-sum_both.t_n(1)]), xlabel('Time [s]', 'FontSize', ftsize), ylab=ylabel('Imp', 'FontSize', ftsize, 'Rotation', 0);
        set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); set(gca, 'FontSize', ftsize, 'YTick', [])
        title('High quality signal segment', 'FontSize', ftsize)
        %set ylim
        rang = range(-1*sum_both.v_n);
        ylim(-1*fliplr([min(-1*sum_both.v_n(:))-0.1*rang, max(-1*sum_both.v_n(:))+0.1*rang]));
        
        % plot template
        subplot(2,3,6), cla, hold on
        time = 0:(length(avtempl)-1); time = time./sum_both.fs;
        for beat_no = 1 : size(ts,1)
            plot(time, -1*ts(beat_no,:), 'color', 0.7*[1 1 1], 'LineWidth', lwidth2)
        end
        plot(time, -1*avtempl, 'r', 'LineWidth', lwidth1), hold off
        set(gca, 'YTick', [])
        xlabel('Time [s]', 'FontSize', ftsize)
        xlim([0, 3]) % Jan Submission
        xlim([0, time(end)])
        %ylab=ylabel('Imp', 'FontSize', ftsize, 'Rotation', 0);
        %set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
        title(['Template (R^2 = ' num2str(R2,2) ')'], 'FontSize', ftsize)
        set(gca, 'FontSize', ftsize)
        %set ylim
        rang = range(ts(:));
        ylim(-1*fliplr([min(ts(:))-0.1*rang, max(ts(:))+0.1*rang]));
    end
        
        
%     annotation('textbox',[0.5, 0.1, 0.1,0.1],'String',{['R2 = ' num2str(R2, 2)] , ['prop breaths bad = ' num2str(prop_bad_breaths, 2) '%'], ['prop dur good = ' num2str(prop_norm_dur,2) '%'], ['norm SD durations = ' num2str(R2min,2) '%']}, 'FontSize', ftsize, 'LineStyle', 'None')
%     if qual
%         annotation('textbox',[0.75, 0.2, 0.1,0.1],'String','High Quality', 'Color', 'b', 'FontSize', ftsize, 'LineStyle', 'None')
%     else
%         annotation('textbox',[0.75, 0.2, 0.1,0.1],'String','Low Quality', 'Color', 'r', 'FontSize', ftsize, 'LineStyle', 'None')
%     end
    savepath = [up.paths.plots_save_folder, save_name];
    %PrintFigs(gcf, paper_size, savepath, up)
    %close all
end


end

function rr_cto = ref_cto(sum_both, up)

% identify peaks
diffs_on_left_of_pt = diff(sum_both.v_n); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt>0);
diffs_on_right_of_pt = diff(sum_both.v_n); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt<0);
peaks = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
% identify troughs
diffs_on_left_of_pt = diff(sum_both.v_n); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt<0);
diffs_on_right_of_pt = diff(sum_both.v_n); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt>0);
troughs = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;
% define threshold
q3 = quantile(sum_both.v_n(peaks), 0.75);
thresh = 0.2*q3;
% find relevant peaks and troughs
extrema = sort([peaks(:); troughs(:)]);
rel_peaks = peaks(sum_both.v_n(peaks) > thresh);
rel_troughs = troughs(sum_both.v_n(troughs) < 0);

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
        cycle_durations(peak_no) = sum_both.t_n(rel_peaks(peak_no+1)) - sum_both.t_n(rel_peaks(peak_no));
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

function rr_wch = ref_wch(sum_both, up)

segLen = 2^nextpow2(12*sum_both.fs);
noverlap = segLen/2;
[data.power, data.freqs] = pwelch(sum_both.v_n,segLen,noverlap, [], sum_both.fs);

% Find spectral peak
[rr_wch, ~, ~] = find_spectral_peak(data, up);

end

function [v, f, p] = find_spectral_peak(data, up)
freq_range = up.paramSet.rr_range/60;

cand_els = zeros(length(data.power),1);
for s = 2 : (length(data.power)-1)
    if data.power(s) > data.power(s-1) & data.power(s) > data.power(s+1) & data.freqs(s) > freq_range(1) & data.freqs(s) < freq_range(2)
        cand_els(s) = 1;
    end
end
clear s freq_range

temp = data.power - min(data.power);  % because the FFT and ACF have -ve values for the power spectra.
[~, r_el] = max(temp.*cand_els); clear cand_els
r_freq = data.freqs(r_el); clear r_el
v = 60*r_freq;

%% Store spectrum
f = data.freqs;
p = data.power;
end

function create_algorithms(up)

fprintf('\n--- Creating logistic regression algorithm ');

%% Load data
loadpath = [up.paths.data_save_folder, up.paths.filenames.data_mat, '.mat'];
load(loadpath);
% use the training dataset to create algorithms
matrices = data_mat.training;

%% Define the response variable
% the manual annotations of high or low quality in the training dataset
resp_var = 'Ann_qual';
rel_col = find(strcmp(matrices.response_header, resp_var));
resp = matrices.response_data(:,rel_col);

%% Define the predictor variables
% don't want to use ID:
rel_cols = find(~strcmp(matrices.design_header, 'id'));
pred = matrices.design_data(:,rel_cols);
pred_names = matrices.design_header(rel_cols);

% Model 2: linear logistic regression (all params)
tbl_names = pred_names; tbl_names{end+1} = resp_var;
ds = mat2dataset([pred,resp],'VarNames',tbl_names);
name = 'lin_log_reg';
mdl = fitglm(ds,'linear','Distribution','binomial','Link','logit');
eval(['algorithms.' name ' = mdl;']);

% prop_norm_dur > 50 && prop_bad_breaths < 10 && R2 >= 0.8

%% Save results matrix to file
save_name = up.paths.filenames.algorithms;
savepath = [up.paths.data_save_folder, up.paths.filenames.algorithms, '.mat'];
save_or_append_data

end

function AnalysePerformances(up)

fprintf('\n--- Analysing Performances')
save_path = [up.paths.data_save_folder, up.paths.filenames.evaluation, '.mat'];
loadpath = [up.paths.data_save_folder, up.paths.filenames.data_mat, '.mat'];
load(loadpath);
data_mat = data_mat.all;

% Check to see whether annotations have been made
rel_col = strcmp(data_mat.response_header, 'rr_ann');
if sum(~isnan(data_mat.response_data(:,rel_col))) == 0
    fprintf('\n - Didn''t analyse performance as no annotations were available')
    return
end

%% Extract required data
vars = {'flat_line', 'rr_novel', 'sqi_novel', 'sqi_agree', 'rr_agree', 'rr_mon', 'win_no', 'subj'};
for var_no = 1 : length(vars)
    curr_var = vars{var_no};
    rel_col = find(strcmp(data_mat.response_header, curr_var));
    eval([curr_var ' = data_mat.response_data(:,rel_col);']);
end
clear curr_var var_no rel_col vars

%% Identify rel windows for each monitoring strategy
rel_windows.novelSQI = sqi_novel == 1;
rel_windows.agreeSQI = sqi_agree == 1;
rel_windows.monitor = ~isnan(rr_mon);

%% Assess proportion of time for which RRs were available
prop_flat_line = 100*mean(flat_line);
no_wins_without_flat_line = sum(flat_line ~=1);
rel_window_types = fieldnames(rel_windows);
subjs = unique(subj);
for rel_window_type_no = 1 : length(rel_window_types)
    curr_rel_window_type = rel_window_types{rel_window_type_no};
    eval(['curr_rel_els = rel_windows.' curr_rel_window_type ';']);
    
    % entire dataset
    eval(['prop_hq.' curr_rel_window_type '.overall = 100*sum(curr_rel_els)/no_wins_without_flat_line;']);
    eval(['prop_hq.' curr_rel_window_type '.overall_n = no_wins_without_flat_line;']);
    
    % per patient basis
    for curr_subj = subjs(:)'
        eval(['curr_subj_rel_els = rel_windows.' curr_rel_window_type ' & subj == curr_subj;']);
        no_subj_wins_without_flat_line = sum(~flat_line & subj == curr_subj);
        temp(curr_subj) = 100*sum(curr_subj_rel_els)/no_subj_wins_without_flat_line;
    end
    eval(['prop_hq.' curr_rel_window_type '.per_subj_lq = quantile(temp, 0.25);']);
    eval(['prop_hq.' curr_rel_window_type '.per_subj_uq = quantile(temp, 0.75);']);
    eval(['prop_hq.' curr_rel_window_type '.per_subj_med = quantile(temp, 0.5);']);
    clear temp no_subj_wins_without_flat_line curr_subj_rel_els curr_subj
    
end
clear rel_window_type_no subjs curr_rel_window_type curr_rel_els rel_windows


%% Identify rel windows for each monitoring strategy
rel_windows.novelSQI_algorithmRR = sqi_novel == 1 & ~isnan(rr_novel);
rel_windows.novelSQI_monitorRR = sqi_novel == 1 & ~isnan(rr_mon);
rel_windows.agreeSQI_algorithmRR = sqi_agree == 1 & ~isnan(rr_novel);
rel_windows.agreeSQI_monitorRR = sqi_agree == 1 & ~isnan(rr_mon);
rel_windows.monitor = ~isnan(rr_mon);

%% Assess proportion of time for which RRs were available
prop_flat_line = 100*mean(flat_line);
no_wins_without_flat_line = sum(flat_line ~=1);
rel_window_types = fieldnames(rel_windows);
subjs = unique(subj);
for rel_window_type_no = 1 : length(rel_window_types)
    curr_rel_window_type = rel_window_types{rel_window_type_no};
    eval(['curr_rel_els = rel_windows.' curr_rel_window_type ';']);
    
    % entire dataset
    eval(['prop_rrs.' curr_rel_window_type '.overall = 100*sum(curr_rel_els)/no_wins_without_flat_line;']);
    eval(['prop_rrs.' curr_rel_window_type '.overall_n = no_wins_without_flat_line;']);
    
    % per patient basis
    for curr_subj = subjs(:)'
        eval(['curr_subj_rel_els = rel_windows.' curr_rel_window_type ' & subj == curr_subj;']);
        no_subj_wins_without_flat_line = sum(~flat_line & subj == curr_subj);
        temp(curr_subj) = 100*sum(curr_subj_rel_els)/no_subj_wins_without_flat_line;
    end
    eval(['prop_rrs.' curr_rel_window_type '.per_subj_lq = quantile(temp, 0.25);']);
    eval(['prop_rrs.' curr_rel_window_type '.per_subj_uq = quantile(temp, 0.75);']);
    eval(['prop_rrs.' curr_rel_window_type '.per_subj_med = quantile(temp, 0.5);']);
    clear temp no_subj_wins_without_flat_line curr_subj_rel_els curr_subj
    
end
clear rel_window_type_no subjs curr_rel_window_type curr_rel_els

%% Assess proportion of time for which an RR was available within the previous five mins
clear rel_els
for rel_window_type_no = 1 : length(rel_window_types)
    
    % Identify relevant windows
    curr_rel_window_type = rel_window_types{rel_window_type_no};
    eval(['curr_rel_els = rel_windows.' curr_rel_window_type ';']);
    
    % Calculate time elapsed since previous RR for each window
    time_elapsed = nan(length(sqi_novel),1);
    prev_subj = 0;
    for s = 1 : length(win_no)
        curr_subj = subj(s);
        curr_win_no = win_no(s);
        if curr_subj ~= prev_subj
            prev_win_no = nan;
        end
        
        % time since previous RR
        if ~flat_line(s)
            time_elapsed(s) = up.paramSet.winLeng*(curr_win_no - prev_win_no);
        end
        
        if curr_rel_els(s) == 1
            prev_win_no = curr_win_no;
        end
        prev_subj = curr_subj;
    end
    clear prev_subj prev_win_no curr_rel_els s curr_subj curr_win_no prev_subj
    
    % Calculate statistics
    rel_els = ~isnan(time_elapsed);
    temp.med = prctile(time_elapsed(rel_els), 50);
    temp.lq = prctile(time_elapsed(rel_els), 25);
    temp.uq = prctile(time_elapsed(rel_els), 75);
    temp.five_min = 100*mean(time_elapsed(rel_els)<5*60);
    temp.ten_min = 100*mean(time_elapsed(rel_els)<10*60);
    
    % store statistics
    eval(['prop_time_elapsed.' curr_rel_window_type ' = temp;']);
    
    clear temp rel_els time_elapsed
end
clear curr_rel_window_type


%% Compare annotated rrs to est rrs (all data deemed to be high qual by each approach - although haven't annotated all monitor windows)

do_section = 1;
if do_section
    
    % identify relevant data to be used for assessing RR performance
    rel_col = find(strcmp(data_mat.response_header, 'for_rr_perf_analysis'));
    rel_rows = data_mat.response_data(:,rel_col) == 1;
    
    % select response data
    rel_col = find(strcmp(data_mat.response_header, 'sqi_novel'));
    novel_sqi = data_mat.response_data(rel_rows,rel_col);
    rel_col = find(strcmp(data_mat.response_header, 'sqi_agree'));
    agree_sqi = data_mat.response_data(rel_rows,rel_col);
    
    % load data
    rel_col = find(strcmp(data_mat.response_header, 'rr_ann'));
    ann_rr = data_mat.response_data(rel_rows,rel_col);
    rel_col = find(strcmp(data_mat.response_header, 'rr_novel'));
    novel_rr = data_mat.response_data(rel_rows,rel_col);
    rel_col = find(strcmp(data_mat.response_header, 'rr_agree'));
    rel_col = find(strcmp(data_mat.response_header, 'rr_novel'));  %%% ADDED IN
    agree_rr = data_mat.response_data(rel_rows,rel_col);
    rel_col = find(strcmp(data_mat.response_header, 'rr_mon'));
    mon_rr = data_mat.response_data(rel_rows,rel_col);
    rel_col = find(strcmp(data_mat.response_header, 'subj'));
    subjs = data_mat.response_data(rel_rows,rel_col);
    subjects = unique(subjs); subjects = subjects(:)';
    clear novel_rrs agree_rrs
    
    % find bias and precision of novel SQI
    counter=0;
    for subj_id = subjects
        rel_els = subjs == subj_id & ~isnan(ann_rr) & ~isnan(novel_rr) & novel_sqi == 1;
        if sum(rel_els) == 0
            continue
        end
        counter = counter+1;
        novel_ids(counter) = unique(subjs(rel_els));
        novel_rrs{counter}(1,:) = ann_rr(rel_els);
        novel_rrs{counter}(2,:) = novel_rr(rel_els);
    end
    [novel_BA, ref_rrs, est_rrs] = perform_BA(novel_rrs);
    savepath = [up.paths.plots_save_folder, 'novel_SQI'];
    make_plots2(ref_rrs, est_rrs, novel_BA, savepath, up);
    novel_BA.no_orig_ids = length(unique(novel_ids));
    
    % find bias and precision of agree SQI
    counter=0;
    for subj_id = subjects
        rel_els = subjs == subj_id & ~isnan(ann_rr) & ~isnan(agree_rr) & agree_sqi == 1;
        if sum(rel_els) == 0
            continue
        end
        counter = counter+1;
        agree_ids(counter) = unique(subjs(rel_els));
        agree_rrs{counter}(1,:) = ann_rr(rel_els);
        agree_rrs{counter}(2,:) = agree_rr(rel_els);
    end
    [agree_BA, ref_rrs, est_rrs] = perform_BA(agree_rrs);
    savepath = [up.paths.plots_save_folder, 'agree_SQI'];
    make_plots2(ref_rrs, est_rrs, agree_BA, savepath, up);
    agree_BA.no_orig_ids = length(unique(agree_ids));
    
    % find bias and precision of monitor
    counter=0;
    for subj_id = subjects
        rel_els = subjs == subj_id & ~isnan(ann_rr) & ~isnan(mon_rr);
        if sum(rel_els) == 0
            continue
        end
        counter = counter+1;
        mon_ids(counter) = unique(subjs(rel_els));
        mon_rrs{counter}(1,:) = ann_rr(rel_els);
        mon_rrs{counter}(2,:) = mon_rr(rel_els);
    end
    [mon_BA, ref_rrs, est_rrs] = perform_BA(mon_rrs);
    savepath = [up.paths.plots_save_folder, 'mon_'];
    make_plots2(ref_rrs, est_rrs, mon_BA, savepath, up);
    mon_BA.no_orig_ids = length(unique(mon_ids));
    clear rel_els subj_id counter ref_rrs est_rrs savepath mon_ids mon_rrs
    
    % find bias and precision of monitor (only in periods when novel SQI showed high quality)
    counter=0;
    for subj_id = subjects
        rel_els = subjs == subj_id & ~isnan(ann_rr) & ~isnan(mon_rr) & novel_sqi == 1;
        if sum(rel_els) == 0
            continue
        end
        counter = counter+1;
        monNov_ids(counter) = unique(subjs(rel_els));
        monNov_rrs{counter}(1,:) = ann_rr(rel_els);
        monNov_rrs{counter}(2,:) = mon_rr(rel_els);
    end
    [monNov_BA, ref_rrs, est_rrs] = perform_BA(monNov_rrs);
    savepath = [up.paths.plots_save_folder, 'novel_SQI_mon'];
    make_plots2(ref_rrs, est_rrs, monNov_BA, savepath, up);
    monNov_BA.no_orig_ids = length(unique(monNov_ids));
    
    %% Create RR results table
    
    Algorithm = {'novel'; 'agree'; 'mon'; 'monNov'};
    for alg_no = 1 : length(Algorithm)
        
        eval(['curr_BA_results = ' Algorithm{alg_no} '_BA;'])
        
        no_wins(alg_no,1) = curr_BA_results.no_samps;
        cp1(alg_no,1) = curr_BA_results.prop_inacc1;
        cp2(alg_no,1) = curr_BA_results.prop_inacc2;
        cp3(alg_no,1) = curr_BA_results.prop_inacc3;
        icp5(alg_no,1) = curr_BA_results.icp5;
        bias_val(alg_no,1) = curr_BA_results.bias.val;
        bias_lci(alg_no,1) = curr_BA_results.bias.lci;
        bias_uci(alg_no,1) = curr_BA_results.bias.uci;
        twoSD_val(alg_no,1) = curr_BA_results.twoSD.val;
        twoSD_lci(alg_no,1) = curr_BA_results.twoSD.lci;
        twoSD_uci(alg_no,1) = curr_BA_results.twoSD.uci;
        prec_val(alg_no,1) = curr_BA_results.prec.val;
        prec_lci(alg_no,1) = curr_BA_results.prec.lci;
        prec_uci(alg_no,1) = curr_BA_results.prec.uci;
        mae(alg_no,1) = curr_BA_results.mae;
        no_orig_ids(alg_no,1) = curr_BA_results.no_orig_ids;
        
    end
    
    rr_results_table = table(Algorithm, no_orig_ids, no_wins, cp1, cp2, cp3, icp5, bias_val, bias_lci, bias_uci, twoSD_val, twoSD_lci, twoSD_uci, prec_val, prec_lci, prec_uci, mae);
    clear agreeNov_* agree_BA agree_ids agree_rr agree_rrs agree_sqi alg_no Algorithm ann_rr bias_lci bias_uci bias_val counter cp1 cp2 cp3 curr_BA_results icp5 mae mon_* no_orig_ids no_wins novel_* orig_id prec_* pred pred_names rel_col rel_els resp sen subj_id subjects subjs monNov_* monNonNov_*
    
    save(save_path, 'rr_results_table');
    clear Algorithm no_wins cp1 cp2 cp3 icp5 bias_val bias_lci bias_uci prec_val prec_lci prec_uci mae
end

%% Output methods:
subjs = unique(data_mat.response_data(:,1));
rel_col = strcmp(data_mat.response_header, 'for_rr_perf_analysis');
for s = 1 : length(subjs), curr_subj = subjs(s); curr_els = data_mat.response_data(:,1) == curr_subj; no_segs(s) = sum(data_mat.response_data(curr_els,rel_col)); end
fprintf('\n - %d annotated segments in total:', sum(no_segs))
fprintf('\n    - %d subjects with 5 segments', sum(no_segs == 5))
fprintf('\n    - %d subjects with 1-5 segments', sum(no_segs >= 1 & no_segs <5))
fprintf('\n    - %d subjects with 0 segments', sum(no_segs == 0))

%% Output results:

fprintf('\n ~~~ RR Performance')
rel_row = find(strcmp(rr_results_table.Algorithm, 'novel'));
fprintf('\n - LOAs of %.1f +/- %.1f bpm when using novel SQI and alg', rr_results_table.bias_val(rel_row), rr_results_table.twoSD_val(rel_row))
rel_row = find(strcmp(rr_results_table.Algorithm, 'monNov'));
fprintf('\n - LOAs of %.1f +/- %.1f bpm when using novel SQI and clinical monitor RRs', rr_results_table.bias_val(rel_row), rr_results_table.twoSD_val(rel_row))
rel_row = find(strcmp(rr_results_table.Algorithm, 'novel'));
fprintf('\n - Frequencies of erroneous RRs: iCP2 %.1f, iCP5 %.1f: when using novel SQI and alg', 100-rr_results_table.cp2(rel_row), rr_results_table.icp5(rel_row))
rel_row = find(strcmp(rr_results_table.Algorithm, 'monNov'));
fprintf('\n - Frequencies of erroneous RRs: iCP2 %.1f, iCP5 %.1f: when using novel SQI and alg', 100-rr_results_table.cp2(rel_row), rr_results_table.icp5(rel_row))
rel_rows = [find(strcmp(rr_results_table.Algorithm, 'novel')), find(strcmp(rr_results_table.Algorithm, 'monNov'))];
rr_results_table(rel_rows,:)
fprintf('\n')

fprintf(' ~~~ Real-World Performance')
fprintf('\n - %.1f %% of the %d non-flat-line segments identified as high quality',prop_hq.novelSQI.overall, prop_hq.novelSQI.overall_n)
fprintf('\n - The %d flat-line segments were excluded from the analysis', sum(flat_line))
fprintf('\n - %.1f (%.1f - %.1f) %% of each subject''s non-flat-line segments identified as high quality',prop_hq.novelSQI.per_subj_med, prop_hq.novelSQI.per_subj_lq, prop_hq.novelSQI.per_subj_uq)
fprintf('\n - RRs estimated in all but %d of the %d high quality segments', prop_hq.novelSQI.overall_n-prop_rrs.novelSQI_algorithmRR.overall_n, prop_hq.novelSQI.overall_n)
fprintf('\n - %d (%d - %d) s between consecutive high quality segments identified by the novel SQI', prop_time_elapsed.novelSQI_algorithmRR.med, prop_time_elapsed.novelSQI_algorithmRR.lq, prop_time_elapsed.novelSQI_algorithmRR.uq)
fprintf('\n - The most recent RR was less than five and ten minutes ago for %.1f and %.1f %% of the recording time', prop_time_elapsed.novelSQI_algorithmRR.five_min, prop_time_elapsed.novelSQI_algorithmRR.ten_min)
fprintf('\n')

% %% Create B-A Summary plots
% 
% % setup
% lwidth = 2;
% ftsize = 16;
% markersize = 7;
% rel_data = rr_results_table.all;
% 
% paper_size = [600,270];
% figure('Position', [20,20,paper_size])
% subplot('Position', [0.25,0.17,0.71,0.81])
% pwv_names = rel_data.Algorithm;
% pwv_axis_labels = strrep(pwv_names, 'novel', 'novel SQI');
% pwv_axis_labels = strrep(pwv_axis_labels, 'agree', 'agree SQI');
% pwv_axis_labels = strrep(pwv_axis_labels, 'mon', 'clinical monitor');
% ba.bias = rel_data.bias_val;
% ba.sd = 0.5*rel_data.prec_val;
% 
% [~,order] = sort(ba.sd, 'ascend');
% 
% counter = 0;
% for pwv_no = order(:)'
%     
%     % Plot the results for this PWV technique
%     counter = counter+1;
%     y_val = counter;
%     line_coords.x_loa = [ba.bias(pwv_no) - 2*ba.sd(pwv_no), ba.bias(pwv_no) + 2*ba.sd(pwv_no)];
%     line_coords.x_bias = ba.bias(pwv_no);
%     line_coords.y = y_val*ones(1,2);
%     % - LOAs
%     plot(line_coords.x_loa, line_coords.y, 'o-k', 'LineWidth', lwidth, 'MarkerFaceColor', 'k', 'MarkerSize', markersize), hold on
%     % - Bias
%     plot(line_coords.x_bias, line_coords.y(1), 'dk', 'LineWidth', lwidth, 'MarkerFaceColor', 'k', 'MarkerSize', markersize+2), hold on
%     
% end
% 
% % Tidy-up
% %xlim([-1,1]*max(abs(xlim)))
% ax = gca;
% %ax.YAxisLocation = 'origin';
% set(gca, 'XGrid','on', 'FontSize', ftsize, 'YTick', 1:length(pwv_axis_labels), 'YTickLabel', pwv_axis_labels)
% xlabel('Error [bpm]', 'FontSize', ftsize)
% ylim([0.5 length(pwv_names)+0.5])
% box off
% legend('95% CI', 'Mean','Location','SouthEast')
% savepath = [up.paths.plots_save_folder, 'BA'];
% PrintFigs(gcf, paper_size/70, savepath)

end

function make_plots(ref_rrs, est_rrs, savepath, up)


%% B-A Plot

% setup
lwidth = 2; ftsize = 22;
paper_size = [9, 6];
figure('Position', [200, 200, 100*paper_size(1), 100*paper_size(2)])
subplot('Position', [0.13, 0.15, 0.85, 0.78])

errors = est_rrs - ref_rrs;
plot(ref_rrs, errors, 'xb', 'LineWidth', lwidth, 'MarkerSize', 8), hold on
xlims = [5*floor(min(ref_rrs)/5), 5*ceil(max(ref_rrs)/5)];
xlims = [5, 40];
plot(xlims,[0,0],'k')
plot(ref_rrs, errors, 'xb', 'LineWidth', lwidth, 'MarkerSize', 8)
xlim(xlims)
xlabel('Reference RR [bpm]', 'FontSize', ftsize)
ylab = ylabel('Error (estimated - reference RR) [bpm]', 'FontSize', ftsize);
set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.47, 0]);
temp.us = strfind(savepath, '_');
temp.slash = strfind(savepath, filesep);
temp.method = savepath(temp.slash(end)+1:temp.us(end)-1);
temp.method = strrep(temp.method, 'novel', 'Novel SQI');
temp.method = strrep(temp.method, 'agree', 'Agreement SQI');
temp.method = strrep(temp.method, 'mon', 'Clinical Monitor');
title(['RR Errors for ' temp.method], 'FontSize', ftsize)
set(gca, 'FontSize', ftsize)
xticks(xlims(1):5:xlims(2))
ylim([-15, 15])
yticks(-20:5:20)
grid on
box off

PrintFigs(gcf, paper_size, [savepath, '_BA'], up)
clear temp

%% Correlation Plot

% setup
lwidth = 2; ftsize = 18;
paper_size = [7, 7];
figure('Position', [200, 200, 100*paper_size(1), 100*paper_size(2)])
subplot('Position', [0.13, 0.15, 0.85, 0.78])

plot(ref_rrs, est_rrs, 'xb', 'LineWidth', lwidth, 'MarkerSize', 8), hold on
temp = [ref_rrs(:); est_rrs(:)];
xlims = [5*floor(min(temp)/5), 5*ceil(max(temp)/5)]; clear temp
xlims = [0, 40];
plot([0,60],[0,60],'k')
plot(ref_rrs, est_rrs, 'xb', 'LineWidth', lwidth, 'MarkerSize', 8)
xlabel('Reference RR [bpm]', 'FontSize', ftsize)
ylab = ylabel('Estimated RR [bpm]', 'FontSize', ftsize);
set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
temp.us = strfind(savepath, '_');
temp.slash = strfind(savepath, filesep);
temp.method = savepath(temp.slash(end)+1:temp.us(end)-1);
temp.method = strrep(temp.method, 'novel', 'Novel SQI');
temp.method = strrep(temp.method, 'agree', 'Agreement SQI');
temp.method = strrep(temp.method, 'mon', 'Clinical Monitor');
title( temp.method, 'FontSize', ftsize)
set(gca, 'FontSize', ftsize)
axis equal
grid on
box off
xlim(xlims)
ylim(xlims)
xticks_to_use = [5*ceil(xlims(1)/5) : 5 : 5*floor(xlims(2)/5)];
xticks(xticks_to_use);
yticks(xticks_to_use);

PrintFigs(gcf, paper_size, [savepath, '_corr'], up)



end

function make_plots2(ref_rrs, est_rrs, novel_BA, savepath, up)

if up.eval.skip_plots == 1
    return
end

% remove segments in which the estimated RRs were zero
rel_els = est_rrs~=0;
ref_rrs = ref_rrs(rel_els);
est_rrs = est_rrs(rel_els);

%% B-A Plot with empirical PDF of errors

% setup
lwidth = 2; ftsize = 22;
paper_size = [9, 6];
figure('Position', [200, 200, 100*paper_size(1), 100*paper_size(2)])
subplot('Position', [0.11, 0.15, 0.61, 0.78])

errors = est_rrs - ref_rrs;
errors(errors<-10) = -10;
errors(errors>10) = 10;
plot(ref_rrs, errors, 'xb', 'LineWidth', lwidth, 'MarkerSize', 8), hold on
xlims = [5*floor(min(ref_rrs)/5), 5*ceil(max(ref_rrs)/5)];
xlims = [5, 42];
% plot(xlims,[0,0],'k')
plot(ref_rrs, errors, 'xb', 'LineWidth', lwidth, 'MarkerSize', 8)
xlim(xlims)
xlabel('Reference RR [bpm]', 'FontSize', ftsize)
ylab = ylabel('Error (estimated - reference RR) [bpm]', 'FontSize', ftsize);
set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.47, 0]);
temp.us = strfind(savepath, '_');
temp.slash = strfind(savepath, filesep);
temp.method = savepath(temp.slash(end)+1:temp.us(end)-1);
temp.method = strrep(temp.method, 'novel_SQI', 'Novel SQI & Clinical Monitor RRs');
temp.method = strrep(temp.method, 'novel', 'Novel SQI & RR Algorithm');
temp.method = strrep(temp.method, 'agree', 'Agreement SQI & RR Algorithm');
temp.method = strrep(temp.method, 'mon', 'Clinical Monitor RRs');
title(temp.method, 'FontSize', ftsize)
set(gca, 'FontSize', ftsize)
xticks(xlims(1):5:xlims(2))
ylims = [-10 10];
if max(errors) > ylims(2)
    error('Change these graph limits');
end
if min(errors) < ylims(1)
    error('Change these graph limits');
end
ylim(ylims)
yticks(-20:5:20)
grid on
box off

% LOAs
plot(xlims,novel_BA.bias.val*ones(1,2),'k')
plot(xlims, novel_BA.lloa.val*ones(1,2), '--k')
plot(xlims, novel_BA.uloa.val*ones(1,2), '--k')
dim = [0.52 0.82 0.1 0.1];
str = {sprintf('Bias:   %.1f', novel_BA.bias.val),sprintf('LOAs: %.1f to %.1f', novel_BA.lloa.val,novel_BA.uloa.val)};
annotation('textbox',dim,'String',str, 'FontSize', ftsize-4);

subplot('Position', [0.79, 0.15, 0.17, 0.78])
histogram(errors, ylims(1):ylims(2))
set(gca,'view',[90 -90])
xlim(ylims)
set(gca, 'FontSize', ftsize, 'XTickLabel', {})
ylab = ylabel('Frequency', 'FontSize', ftsize);
grid on
box off
ylims2 = ylim;
hold on, plot([0,0],ylims2,'k')

PrintFigs(gcf, paper_size, [savepath, '_BA'], up)
clear temp

%% Correlation Plot

% setup
lwidth = 2; ftsize = 18;
paper_size = [7, 7];
figure('Position', [200, 200, 100*paper_size(1), 100*paper_size(2)])
subplot('Position', [0.13, 0.15, 0.85, 0.78])

plot(ref_rrs, est_rrs, 'xb', 'LineWidth', lwidth, 'MarkerSize', 8), hold on
temp = [ref_rrs(:); est_rrs(:)];
xlims = [5*floor(min(temp)/5), 5*ceil(max(temp)/5)]; clear temp
xlims = [0, 42];
plot([0,60],[0,60],'k')
plot(ref_rrs, est_rrs, 'xb', 'LineWidth', lwidth, 'MarkerSize', 8)
xlabel('Reference RR [bpm]', 'FontSize', ftsize)
ylab = ylabel('Estimated RR [bpm]', 'FontSize', ftsize);
set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
temp.us = strfind(savepath, '_');
temp.slash = strfind(savepath, filesep);
temp.method = savepath(temp.slash(end)+1:temp.us(end)-1);
temp.method = strrep(temp.method, 'novel_SQI', 'Novel SQI & Clinical Monitor RRs');
temp.method = strrep(temp.method, 'novel', 'Novel SQI & RR Algorithm');
temp.method = strrep(temp.method, 'agree', 'Agreement SQI & RR Algorithm');
temp.method = strrep(temp.method, 'mon', 'Clinical Monitor RRs');
title(temp.method, 'FontSize', ftsize)
set(gca, 'FontSize', ftsize)
axis equal
grid on
box off
xlim(xlims)
ylim(xlims)
xticks_to_use = [5*ceil(xlims(1)/5) : 5 : 5*floor(xlims(2)/5)];
xticks(xticks_to_use);
yticks(xticks_to_use);

PrintFigs(gcf, paper_size, [savepath, '_corr'], up)



end

function PrintFigs(h, paper_size, savepath, up, do_close)
set(h,'PaperUnits','inches');
set(h,'PaperSize', [paper_size(1), paper_size(2)]);
set(h,'PaperPosition',[0 0 paper_size(1) paper_size(2)]);
print(h,'-dpdf',savepath)
print(h,'-depsc',savepath)
%print(h,'-dpng',savepath)

if nargin == 5 && do_close == 0
else
    close all
end
end

function [BA, ref_rrs, est_rrs] = perform_BA(rrs)

% calculate errors (in bpm)
[subjs, errors, ref_rrs, est_rrs] = deal([]);
for subj_no = 1 : length(rrs)
    rel_ref_rrs = rrs{subj_no}(1,:);
    rel_est_rrs =  rrs{subj_no}(2,:);
    
    % remove segments in which estimated RRs were zero
    rel_els = rel_est_rrs~=0;
    rel_ref_rrs = rel_ref_rrs(rel_els);
    rel_est_rrs = rel_est_rrs(rel_els);
    
    rel_errors = rel_ref_rrs - rel_est_rrs;
    errors = [errors, rel_errors];
    ref_rrs = [ref_rrs, rel_ref_rrs];
    est_rrs = [est_rrs, rel_est_rrs];
    subjs = [subjs, subj_no*ones(1,length(rel_errors))];
end

% calculate MAE
BA.mae = nanmean(abs(errors));

% calculate bias
BA.bias.val = nanmean(errors);
BA.no_samps = sum(~isnan(errors));
BA.prop_inacc1 = 100*sum(~isnan(errors) & abs(errors) < 1)/sum(~isnan(errors));
BA.prop_inacc2 = 100*sum(~isnan(errors) & abs(errors) < 2)/sum(~isnan(errors));
BA.prop_inacc3 = 100*sum(~isnan(errors) & abs(errors) < 3)/sum(~isnan(errors));
BA.cp2 = 100*sum(~isnan(errors) & abs(errors) < 2)/sum(~isnan(errors));
BA.icp5 = 100*sum(~isnan(errors) & abs(errors) > 5)/sum(~isnan(errors));

% calculate bias CI
n = length(rrs);            % sample size (no subjects)
doff = n - 1;               % Degrees of Freedom
tval = tinv(0.975,doff);    % value from T distribution
std_diffs = std(errors(~isnan(errors)));    % standard dev of diffs
std_err_diffs = sqrt((std_diffs^2)/n);      % standard err of diffs
BA.bias.lci = BA.bias.val - (tval*std_err_diffs);
BA.bias.uci = BA.bias.val + (tval*std_err_diffs);

% calculate loa
[p,tbl] = anova1(errors,subjs,'off'); % one-way anova of the differences
MS_subj = tbl{2,4};
MS_resid = tbl{3,4};
no_wins = [];
for subj_no = 1 : length(rrs)
    no_wins(subj_no) = sum(subjs == subj_no);
end
denominator = ( (sum(no_wins)^2) - sum(no_wins.^2) ) / ( (n-1)*sum(no_wins) );
heterogeneity_variance = (MS_subj - MS_resid)/denominator;
total_variance = heterogeneity_variance + MS_resid;
BA.lloa.val = BA.bias.val - (1.96*sqrt(total_variance));
BA.uloa.val = BA.bias.val + (1.96*sqrt(total_variance));

% calculate loa CI
std_err_lim = sqrt((3*(std_diffs^2))/n);      % standard err of the limit of diffs
BA.lloa.lci = BA.lloa.val - (tval*std_err_lim);
BA.lloa.uci = BA.lloa.val + (tval*std_err_lim);
BA.uloa.lci = BA.uloa.val - (tval*std_err_lim);
BA.uloa.uci = BA.uloa.val + (tval*std_err_lim);

% calculate 2SD
BA.twoSD.val = 1.96*sqrt(total_variance);
BA.twoSD.lci = BA.twoSD.val - (tval*std_err_lim);
BA.twoSD.uci = BA.twoSD.val + (tval*std_err_lim);

% calculate precision
BA.prec.val = BA.uloa.val - BA.lloa.val;
BA.prec.lci = BA.prec.val - (2*std_err_lim*tval);
BA.prec.uci = BA.prec.val + (2*std_err_lim*tval);

end