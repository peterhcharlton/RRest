function [overall_peaks, overall_onsets] = co_ppg_peak_detector(ppg, fs, up)
%CO_PPG_PEAK_DETECTOR detects PPG pulse peaks
%
% This has been adapted by Peter Charlton from code supplied by Christina
% Orphanidou.
%
% Purpose: To detect PPG peaks.
%
% With grateful thanks to Christina Orphanidou and Alexander Darrell.
%

%% Pre-processing
% remove nans
ppg.v(isnan(ppg.v)) = mean(ppg.v(~isnan(ppg.v)));
% filter
ppgfilt.v = filtfilt(up.paramSet.CO_peak_det.PPG_smoother.coeffs, 1, ppg.v);
ppgfilt.t = ppg.t;

%% Setup
% Segmentation into 10s windows
win_length = 10;
overlap = 3;
win_starts = ppgfilt.t(1):win_length:ppgfilt.t(end);
win_ends = win_starts + win_length + overlap;
win_ends(end) = win_ends(end) - overlap;
% setup variables
[overall_peaks, overall_onsets] = deal([]);

%% Perform peak detection on each 10s window
for win_no = 1 : length(win_starts)
    
    rel_data = struct;
    
    rel_data.i = find(ppgfilt.t >= win_starts(win_no) & ppgfilt.t <= win_ends(win_no));
    rel_data.t = ppgfilt.t(rel_data.i);
    rel_data.v = ppgfilt.v(rel_data.i);
    
    %% Calculate thresholds
    thresh1=quantile(rel_data.v,up.paramSet.CO_peak_det.upctl);
    thresh2=quantile(rel_data.v,up.paramSet.CO_peak_det.lpctl);
    thresh3=thresh2+0.3*(thresh1-thresh2);
    thresh4=thresh2+0.7*(thresh1-thresh2);
    
    %% Find all peaks
    % identify peaks
    diffs_on_left_of_pt = diff(rel_data.v); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt>0);
    diffs_on_right_of_pt = diff(rel_data.v); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); level_on_right_of_pt = logical(diffs_on_right_of_pt==0); diffs_on_right_of_pt = logical(diffs_on_right_of_pt<0);
    diffs_on_second_right_of_pt = [1; diff(rel_data.v(2:end))]; diffs_on_second_right_of_pt = diffs_on_second_right_of_pt(2:end); diffs_on_second_right_of_pt = logical(diffs_on_second_right_of_pt<0);
    peaks.i = find((diffs_on_left_of_pt & diffs_on_right_of_pt) | ...
        (diffs_on_left_of_pt & level_on_right_of_pt & diffs_on_second_right_of_pt))+1;
    peaks.i = peaks.i';
    % Take maximum value in window if no peaks were found
    if isempty(peaks.i)
        peaks.i = find(max(rel_data.v) == rel_data.v);
    end
    % Extract signal values at peaks
    peaks.t = rel_data.t(peaks.i);
    peaks.v = rel_data.v(peaks.i);
    
    %% Classify peaks
    % Identify relevant peaks according to amplitude
    upperdiff = abs(peaks.v-thresh1);
    middlehighdiff = abs(peaks.v-thresh4);
    middlelowdiff = abs(peaks.v-thresh3);
    lowerdiff = abs(peaks.v-thresh2);
    upper_pks = find(upperdiff<middlehighdiff & upperdiff<middlelowdiff & upperdiff<lowerdiff);
    PPG_PKS.i = peaks.i(upper_pks);
    PPG_PKS.v = peaks.v(upper_pks);
    PPG_PKS.t = peaks.t(upper_pks);
    % eliminate peaks which are too close together in time
    good_els = find(diff(PPG_PKS.i)>=fs/3)+1;
    PPG_PKS.i = PPG_PKS.i(good_els);
    PPG_PKS.v = PPG_PKS.v(good_els);
    PPG_PKS.t = PPG_PKS.t(good_els);
    
    %% Find troughs
    
    PPG_TRS.i = nan(length(PPG_PKS.t)-1,1);
    for s = 1 : (length(PPG_PKS.t)-1)
        start_el = PPG_PKS.i(s);
        [~, additional_el] = min(rel_data.v(PPG_PKS.i(s):PPG_PKS.i(s+1)));
        PPG_TRS.i(s) = start_el - 1 + additional_el;
    end
    PPG_TRS.v = rel_data.v(PPG_TRS.i);
    PPG_TRS.t = rel_data.t(PPG_TRS.i);
    
    temp_pk_is = rel_data.i(PPG_PKS.i(:));
    temp_on_is = rel_data.i(PPG_TRS.i(:));
    
    overall_peaks = [overall_peaks; temp_pk_is(:)]; overall_peaks = unique(overall_peaks);
    overall_onsets = [overall_onsets; temp_on_is(:)]; overall_onsets = unique(overall_onsets);
    
    clear upper lower middlehigh middlelow temp_pk_is temp_on_is PPG_TRS PPG_PKS additional_el s good_els peaks lower_pks upper_pks lowerdiff middlelowdiff middlehighdiff upperdiff ind
    
end
