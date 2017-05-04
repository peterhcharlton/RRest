function calculate_sqi(up)
%CALCULATE_SQI calculates SQIs from the filtered ECG and PPG signals
%	            calculate_sqi(up)
%
%   Reference:
%       Orphanidou, C. et al., 2015.
%       Signal-quality indices for the electrocardiogram and photoplethysmogram: derivation and applications to wireless monitoring.
%       IEEE Journal of Biomedical and Health Informatics, 19(3), pp.832–8.
%       Available at: http://www.ncbi.nlm.nih.gov/pubmed/25069129.
%
%	Inputs:
%		data            data files should be stored in the specified format
%       up              universal parameters structure
%
%	Outputs:
%       	for each subject, n:
%       n_sqi.m         - a file of SQI values
%
% With grateful thanks to Christina Orphanidou for supplying some of the
% code contained within this stage.

fprintf('\n--- Calculating SQIs ');

% For each subject
for subj = up.paramSet.subj_list
    
    %% Skip if this processing has been done previously
    save_name = 'sqi';
    savepath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.sqi, '.mat'];
    exist_log = check_exists(savepath, save_name);
    if exist_log
        continue
    end
    
    %% Load window timings
    loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.win_timings, '.mat'];
    load(loadpath);
    
    %% Cycle through each signal
    for curr_sig_type = {'ekg', 'ppg'}
        sig_type = curr_sig_type{1,1};
        % Cycle through each signal of this type
        eval(['sigs = up.paramSet.' sig_type '_sigs;']);
        for sig_no = 1 : length(sigs)
            curr_sig = sigs{sig_no};
            
            %% Load PPG and ECG signals
            loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.int_respSigs];
            rel_name1 = [curr_sig, up.paths.filenames.elim_vhf];
            load(loadpath, rel_name1)
            eval(['sig_data = ' rel_name1 ';']); clear rel_name1
            %% Load beat detections if CO's peak detector has already been used
            if strcmp(sig_type, 'ppg')
                rel_name2 = [curr_sig, up.paths.filenames.pulse_peaks, 'COr', up.paths.filenames.elim_vhf];
                threshold = 0.86;                                               % ppg cross-correlation threshold value (for sqi)
            else
                rel_name2 = [curr_sig, up.paths.filenames.qrss, up.al.options.RDt{1,1}, up.paths.filenames.elim_vhf];
                threshold = 0.66;                                               % ecg cross-correlation threshold value (for sqi)
            end
            file_contents = whos('-file',loadpath);
            var_names = extractfield(file_contents, 'name');
            if sum(strcmp(var_names, rel_name2))
                load(loadpath, rel_name2);
                eval(['beats = ' rel_name2 ';']); clear rel_name2
            else
                %% Perform PPG beat detection if they don't exist
                % High-pass filter data for peak detection
                s_filt = elim_sub_cardiac(sig_data, up);
                
                % Detect Pulse Peaks in HPF'd version
                [peaks, onsets] = co_ppg_peak_detector(s_filt, s_filt.fs, up);
                beats.p.t = s_filt.t(peaks);
                beats.p.v = s_filt.v(peaks);
                
            end
            
            %% Filter PPG to eliminate sub-cardiac freqs
            if strcmp(sig_type, 'ppg')
                if ~exist('s_filt', 'var')
                    s_filt = elim_sub_cardiac(sig_data, up);
                end
                sig_data.v = s_filt.v; clear s_filt
            end
            
            %% further beat detection information
            beats.rr.t = beats.p.t(2:end);
            beats.rr.v = 1000*diff(beats.p.t);   % beat-to-beat intervals in ms
            beats.p.i = nan(1,length(beats.p.t));
            for s = 1 : length(beats.p.t)
                beats.p.i(s) = find(sig_data.t == beats.p.t(s));
            end
            beats.rr.i = nan(1,length(beats.rr.t));
            for s = 1 : length(beats.rr.t)
                beats.rr.i(s) = find(sig_data.t == beats.rr.t(s));
            end
            clear s
            
            % Setup variables
            [SQI, coeff, hr] = deal(ones(length(wins.t_start),1));
            
            for win_no = 1 : length(wins.t_start)
                % identify 10 s sections within this window:
                win_start = wins.t_start(win_no);
                win_end = wins.t_end(win_no);
                win_length = win_end - win_start;
                
                sec_starts = win_start : 10 : (win_end-10);
                if rem(win_length, 10) > 0
                    sec_starts = [sec_starts, (win_end-10)];
                end
                sec_ends = sec_starts + 10;
                
                [temp_SQI, temp_coeff] = deal(nan(length(sec_starts),1));
                for seg_no = 1 : length(sec_starts)
                    
                    % find relevant data for this section:
                    section_data.i = find(sig_data.t>=sec_starts(seg_no) & ...
                        sig_data.t<=sec_ends(seg_no));
                    section_data.t = sig_data.t(section_data.i);
                    section_data.v = sig_data.v(section_data.i);
                    
                    rel_beats = find(beats.p.t>=sec_starts(seg_no) & ...
                        beats.p.t<sec_ends(seg_no));
                    dummy_beat_inds = beats.p.i(rel_beats);
                    rel_rrs = find(beats.rr.t>=sec_starts(seg_no) & ...
                        beats.rr.t<sec_ends(seg_no));
                    dummy_rr_inds = beats.rr.i(rel_rrs);
                    % If there's one or no beats, put SQI = 0:
                    if isempty(dummy_beat_inds) || isempty(dummy_rr_inds)
                        temp_coeff(seg_no) = nan;
                        temp_SQI(seg_no) = 0;
                        clear dummy_beat_inds section_data rel_beats section_beats section_rr section_length
                        continue
                    end
                    section_beats.v = beats.p.v(rel_beats);
                    section_beats.i = dummy_beat_inds - section_data.i(1);
                    section_beats.t = sig_data.t(dummy_beat_inds);
                    section_rr.i = beats.rr.t(rel_rrs);
                    section_rr.v = beats.rr.v(rel_rrs);
                    
                    % calculate this section's SQI
                    section_length = ceil(section_data.t(end) - section_data.t(1));
                    [temp_coeff(seg_no),temp_SQI(seg_no)]= sqi_calculator(section_data.v, section_beats.i, sig_data.fs, threshold, section_length, up);
                    
                    clear dummy_beat_inds section_data rel_beats section_beats section_rr section_length
                    
                end
                
                if sum(temp_SQI) == length(temp_SQI)
                    SQI(win_no) = 1;
                else
                    SQI(win_no) = 0;
                end
                coeff(win_no) = mean(temp_coeff);
                
                clear  seg_no temp_coeff temp_SQI win_start win_end win_length sec_ends sec_starts
                
                %% Find HR for this section
                rel_beats = find(beats.p.t >= wins.t_start(win_no) & beats.p.t < wins.t_end(win_no));
                if ~isempty(rel_beats)
                    duration_of_beats = beats.p.t(rel_beats(end)) - beats.p.t(rel_beats(1));
                    ave_duration_of_beat = duration_of_beats/(length(rel_beats)-1);
                    hr(win_no) = 60/ave_duration_of_beat;
                else
                    hr(win_no) = nan;
                end
                
            end
            
            % put results into variable to store:
            temp.t = wins.t_start;
            temp.v = SQI;
            temp.hr = hr;
            
            eval(['sqi.' curr_sig ' = temp;']);  clear temp coeff HR SQI sig_data beats threshold win_no
            
        end
        
    end
    
    %% Save SQIs to file
    save_or_append_data
    clear ekg* ppg* wins sig sqi
    
end

end

function [R2,SQI] = sqi_calculator(data,beats,fs,IN,w, up)
%calculates SQI for single segment of data.
% taken from CO's code. Adapted by PC 21/06/2013
%Inputs:    data - 1D timeseries of data (ecg or ppg)
%           beats - indices of beats (qrs spikes or ppg pulses)
%           fs - sampling rate
%           IN - threshold for average correlation coefficient (default 0.66 ecg, 0.86 ppg)
%           w - window size (default: 10s).
%Outputs:   R2 - average correlation coefficient
%           SQI - 0 if bad and 1 if good

% set up variables:
SQI=0; R2=0; avtempl=0; ts=0;

%find mean RR interval to define size of template
hrs=floor(mean(diff(beats)));
q=[0 beats fs*10];

%% Apply rules: Rule 1 (40<HR) || Rule 1 (HR<180) || Rule 2 (are there any gaps > 3s?) || Rule 3 (ratio of max to min RR interval should be less than 2.2)
if length(beats) < (40*w/60) || length(beats) > (180*w/60) || ...
        isempty(find(diff(q)>3*fs))==0 || max(diff(peaks))/min(diff(peaks))>2.2
    return
else
    %% template matching-first identify QRS complexes which have a complete template within the window of data.
    hr=60*fs./hrs;
    ts=[];
    j=find(beats>hrs/2);
    l=find(beats+floor(hrs/2)<length(data));
    if isempty(l)==1
        return
    else
        %find QRS complexes
        for k=j(1):l(end)
            t=data(beats(k)-floor(hrs/2):beats(k)+floor(hrs/2));
            tt=t/norm(t); tt = tt(:)';
            ts=[ts;tt];
        end
    end
end

%find all templates in current window
avtempl=mean(ts,1);

%now calculate correlation for every beat in this window
r2 = nan(size(ts,1),1);
for k=1:size(ts,1)
    r2(k)=corr2(avtempl,ts(k,:));
end

%calculate mean correlation coefficient
R2=mean(abs(r2));

if R2<IN
    SQI=0;
else
    SQI=1;
end

%% Plot template and inidividual beats
plot_ex = 0;
if plot_ex
    paper_size = [6, 5];
    figure('Position', [200, 200, 100*paper_size(1), 100*paper_size(2)], 'Color',[1 1 1])
    lwidth1 = 3; lwidth2 = 2; ftsize = 22;
    time = 0:(length(avtempl)-1); time = time./fs;
    hold on,
    for beat_no = 1 : size(ts,1)
        plot(time, ts(beat_no,:), 'color', 0.7*[1 1 1], 'LineWidth', lwidth2)
    end
    plot(time, avtempl, 'r', 'LineWidth', lwidth1)
    set(gca, 'YTick', [])
    xlabel('Time [s]', 'FontSize', ftsize)
    xlim([0, time(end)])
    if IN==0.86
        ylab=ylabel('PPG', 'FontSize', ftsize, 'Rotation', 0);
    else
        ylab=ylabel('ECG', 'FontSize', ftsize, 'Rotation', 0);
    end
    set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    set(gca, 'FontSize', ftsize, 'XTick', [])
    %set ylim
    rang = range(ts(:));
    ylim([min(ts(:))-0.1*rang, max(ts(:))+0.1*rang]);
    if SQI && IN==0.86
        save_name = 'ppg_high_qual';
        annotation('textbox',[0.01, 0.9, 0.1,0.1],'String','(b)', 'FontSize', ftsize, 'LineStyle', 'None')
    elseif SQI && IN==0.66
        save_name = 'ecg_high_qual';
        annotation('textbox',[0.01, 0.9, 0.1,0.1],'String','(a)', 'FontSize', ftsize, 'LineStyle', 'None')
    elseif ~SQI && IN==0.66
        save_name = 'ecg_low_qual';
        annotation('textbox',[0.01, 0.9, 0.1,0.1],'String','(c)', 'FontSize', ftsize, 'LineStyle', 'None')
    elseif ~SQI && IN==0.86
        save_name = 'ppg_low_qual';
        annotation('textbox',[0.01, 0.9, 0.1,0.1],'String','(d)', 'FontSize', ftsize, 'LineStyle', 'None')
    end
    savepath = [up.paths.plots_save_folder, save_name];
    %PrintFigs(gcf, paper_size, savepath, up)
    if ~SQI
        a=1;
    end
    close all
end

end

function PrintFigs(h, paper_size, savepath, up)
set(h,'PaperUnits','inches');
set(h,'PaperSize', [paper_size(1), paper_size(2)]);
set(h,'PaperPosition',[0 0 paper_size(1) paper_size(2)]);
print(h,'-dpdf',savepath)
print(h,'-dpng',savepath)

% you need to download 'export_fig' from:
% http://uk.mathworks.com/matlabcentral/fileexchange/23629-export-fig
export_fig_dir_path = 'C:\Users\pc13\Documents\GitHub\phd\Tools\Other Scripts\export_fig\';
addpath(export_fig_dir_path)
export_fig(savepath, '-eps')

close all
end