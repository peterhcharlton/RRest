function run_mimic_imp_annotation
% run_mimic_imp_annotation is designed to annotate the quality and breaths
% of signals in the impedance pneumography signal quality index (SQI) dataset
% (extracted from the MIMIC III database). It was used to evaluate the
% impedance pneumography SQI (using only the breath annotation functionality).
%
%               run_mimic_imp_annotation
%
%	Inputs:
%       mimic_imp_sqi_data.mat - a file containing curated data from the MIMIC database.
%        - Specify the path to this file in the "setup_universal_params" function below.
%
%	Outputs:
%       Matlab files containing the annotations for each subject, called
%          "ImP_SQI#_##_an.mat" (where #_## is the subject number, followed
%          by the initials of the annotator)
%           
%   Further Information:
%       This version of the run_mimic_imp_annotation is provided to facilitate
%       replication of the analysis reported in:
%           Charlton P. H. et al., "An impedance pneumography signal
%           quality index for respiratory rate monitoring: design,
%           assessment and application", [under review]
%       Further information on this study can be obtained at:
%           http://peterhcharlton.github.io/RRest/imp_sqi.html
%       In addition, further information on RRest, a toolbox of respiratory
%       algorithms which can be used with this dataset, can be obtained at:
%           http://peterhcharlton.github.io/RRest/index.html
%
%   Comments, Questions, Criticisms, Feedback, Contributions:
%       See: http://peterhcharlton.github.io/RRest/contributions.html
%
%   Version:
%       v.0.1 - accompanying peer review, 6th Aug 2020 by Peter H Charlton
%
%   Source:
%       This script has been adapted from 'MIMICII_data_importer.m', which
%       is one of the scripts in the RRest toolbox, available at:
%           https://github.com/peterhcharlton/RRest
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
% ~~~ This function needs adjusting by the individual user as it contains paths
% specific to their computer ~~~
up = setup_universal_params;

%% Initialisation

% Load data
load(up.paths.data_file)
subjects = 1:length(data);   % create a list of subjects which are to be analysed.

% mode
up.mode = 2;        % select 1 if annotating high or low quality, 2 if annotating breaths
if up.mode == 1, fprintf('\n *** Mode selected: quality annotation *** '), elseif up.mode == 2, fprintf('\n *** Mode selected: Breath annotation *** '), end

fprintf(['\n\n\n\n\n~~~~~~~  Starting annotations  ~~~~~~~\n' ...
    'Instructions: Use the left mouse key to\n' ...
    'annotate a peak, and use the right mouse\n' ...
    'key to annotate a trough. To delete the\n' ...
    'most recent annotation, hold down SHIFT\n' ...
    'and press a mouse key.\n' ...
    'Please note that any previous annotations\n' ...
    'will be kept as long as the results files\n' ...
    'are in the correct location.\n\n']);

operator=input('(To enter review mode, type  review  )\n\nOtherwise, please type your first and last initial,\nand press enter:\n','s');
%operator = 'PC';
if strcmp(operator, 'review') == 1
    paw_verification_script_checker(subjects, data, up);
    return
end

question=input('\nWould you like to start from the \nfirst subject? (y or n, followed by enter)\n', 's');
%question = 'y';
if strcmp(question, 'n') == 1
    start_subject = str2num(input('\nWhich subject?\n', 's'));
else
    start_subject = 1;
end

%% Cycle through different subjects and periods
start_subject_el = find(subjects == start_subject);
% load data on which windows were high quality if annotating breaths
if up.mode == 2
    load(up.paths.data_file2)
    rel_col = find(strcmp(data_mat.all.response_header, 'subj'));
    data_mat.all.subj = data_mat.all.response_data(:,rel_col);
    rel_col = find(strcmp(data_mat.all.response_header, 'win_no'));
    data_mat.all.win_no = data_mat.all.response_data(:,rel_col);
end
if isempty(start_subject_el)
    fprintf('\nThis subject doesn''t exist. Try all over again.');
else
    subjects = subjects(start_subject_el:end);
    for s = subjects(1:end)
        
        %% Identify wins deemed to be high qual
        if up.mode == 2 % i.e. identifying RRs
            rel_col1 = find(strcmp(data_mat.all.response_header, 'sqi_novel'));
            rel_els1 = data_mat.all.subj==s;
            rel_els2 = data_mat.all.response_data(:, rel_col1) == 1;
            high_qual_wins = data_mat.all.win_no(rel_els1 & rel_els2);
            if length(high_qual_wins) > 5
                high_qual_wins = high_qual_wins(1:5);
            end
        end
                
        %% Load Impedance
        imp.fs = data(s).ref.resp_sig.imp.fs;
        imp.v = data(s).ref.resp_sig.imp.v;
        imp.t = [0:(length(imp.v)-1)]/imp.fs;
        
        %% Process Impedance
        imp = lpf_to_exclude_resp(imp, up);
        imp.v = -1*imp.v;
        
        %% Set up file for keeping annotations in:
        % see if this subject has already been annotated:
        if up.mode == 1
            qual.t = []; qual.v = [];
            save(up.paths.annotations_file, 'qual', 'up')
        else
            pk_anns.t=[]; pk_anns.v=[];
            save(up.paths.annotations_file, 'pk_anns', 'up')
        end 
        
        %% plot waveforms and annotations (for further annotation)
        duration_of_signal = imp.t(end) - imp.t(1);
        NUMWINS = floor(duration_of_signal / up.win_length); clear duration_of_signal
        edge_time = up.edge_prop*up.win_length;
        for win = 1:NUMWINS
            % skip if annotating breaths and this window was deemed to be
            % low quality
            if up.mode == 2 && ~sum(win==high_qual_wins)
                continue
            end
            
            % Find window timings
            starttime = imp.t(1)+(win-1)*up.win_length;
            endtime = starttime + up.win_length;
            
            % Skip if this window has already been annotated
            savename = [up.paths.annotations_folder, 'ImP_SQI' num2str(s), '_' operator '_an.mat'];
            if exist(savename, 'file')
                temp = load(savename);
                rel_els = temp.pk_anns.t>=starttime & temp.pk_anns.t<=endtime;
                min_diff = min(diff(temp.pk_anns.t(rel_els)));
                if min_diff < 1
                    error(['Check window ', num2str(win), ', subj ', num2str(s)])
                end
                if sum(rel_els)>0
                    continue
                end
            end
            
            % Setup plot
            rel_els = find(imp.t >= starttime & imp.t < endtime);
            grey_els = find(imp.t >= (starttime - edge_time) & imp.t < (endtime + edge_time));
            scrsz = get(0,'ScreenSize');
            load(up.paths.annotations_file)
            
            % Plot
            h=figure('OuterPosition',[1 1 scrsz(3) scrsz(4)]);
            %save('plot_vars.mat', 'h', 'imp', 'Re', 'grey_paw_els', 'rel_paw_els', 'rel_Re_els', 'grey_Re_els', 'period', 's', 'subjects')
            [axis_h, h2] = annotation_plotter2(up.paths.annotations_file, h, imp, grey_els, rel_els, s, subjects, up, win, NUMWINS,starttime);
            
            % use pointer to identify breaths
            set(gcf, 'Pointer', 'cross');
            if up.mode == 1
                % Construct a questdlg with three options
                choice = MFquestdlg([0.1,0.1],'High quality?', ...
                    '', ...
                    'Yes','No','No');
                % Handle response
                switch choice
                    case 'Yes'
                        load(up.paths.annotations_file);
                        [qual.t,inds]=unique([qual.t; starttime]);
                        qual.v=[qual.v; 1]; qual.v = qual.v(inds);
                        [qual.t, inds] = sort(qual.t);
                        qual.v = qual.v(inds);
                        save(up.paths.annotations_file, 'qual');
                    case 'No'
                        load(up.paths.annotations_file);
                        [qual.t,inds]=unique([qual.t; starttime]);
                        qual.v=[qual.v; 0]; qual.v = qual.v(inds);
                        [qual.t, inds] = sort(qual.t);
                        qual.v = qual.v(inds);
                        save(up.paths.annotations_file, 'qual');
                end
            else
                set(h2,'ButtonDownFcn',{@Click_CallBack2 gca axis_h h2 up up.paths.annotations_file starttime});
                pause
            end
            close(h);
            
        end
        
        if up.DOSAVE
            fprintf(['~~~~~~ Saving subject ' num2str(s) '  ~~~~~~\n']);
            if exist(savename, 'file')
                old_data = load(savename);
            end
            new_data.up = up;
            temp = load(up.paths.annotations_file);
            % If there are annotations in both the old file and the current annotations file
            if exist('old_data', 'var') && sum(strcmp(fieldnames(old_data), 'pk_anns')) && ~isempty(old_data.pk_anns.t) ...
                    && sum(strcmp(fieldnames(temp), 'pk_anns')) && ~isempty(temp.pk_anns.t)
                all_t = [temp.pk_anns.t; old_data.pk_anns.t];
                all_v = [temp.pk_anns.v; old_data.pk_anns.v];
                [~, order] = sort(all_t);
                new_data.pk_anns.t = all_t(order);
                new_data.pk_anns.v = all_v(order);
                clear order all_t all_v
            % If there are annotations in the old file
            elseif exist('old_data', 'var') && sum(strcmp(fieldnames(old_data), 'pk_anns')) && ~isempty(old_data.pk_anns.t)
                new_data.pk_anns = old_data.pk_anns;
            % If there are annotations in the current annotations file
            elseif sum(strcmp(fieldnames(temp), 'pk_anns')) && ~isempty(temp.pk_anns.t)
                new_data.pk_anns = temp.pk_anns;
            else
                [new_data.pk_anns.t, new_data.pk_anns.v] =deal([]);
            end
            if exist('qual', 'var') && ~isempty(qual.t)
                new_data.qual = qual;
            elseif sum(strcmp(fieldnames(temp), 'qual'))  && ~isempty(temp.qual.t)
                new_data.qual = temp.qual;
            elseif exist('old_data', 'var') && sum(strcmp(fieldnames(old_data), 'qual'))  && ~isempty(old_data.qual.t)
                new_data.qual = old_data.qual;
            else
                new_data.qual.t = []; qual.v = [];
            end
            pk_anns = new_data.pk_anns;
            qual = new_data.qual;
            up = new_data.up;
            save(savename, 'pk_anns', 'qual', 'up')
            clear old_data new_data temp pk_anns qual
        end
        
        clear pk_anns qual NUMWINS pawmean STUDYID
        
    end
    
end

end

function up = setup_universal_params

close all

%% folderpaths
up.paths.root_data_folder = '/Users/petercharlton/Documents/Data/ImP_SQI_mimic_dataset_test/';
up.paths.data_file = [up.paths.root_data_folder, 'mimic_imp_sqi_data.mat'];
up.paths.data_file2 = [up.paths.root_data_folder, 'data_mat.mat'];
up.paths.resultspath = [up.paths.root_data_folder, 'SQI_development', filesep];
up.paths.annotations_folder = [up.paths.root_data_folder, '2019_annotations', filesep];
up.paths.annotations_file = [up.paths.root_data_folder, '2019_annotations', filesep, 'temp.mat'];
up.paths.dlg = '/Users/petercharlton/Google Drive/Work/Projects/PhD/Github/phd/Tools/Other Scripts/MFquestdlg/';
addpath(up.paths.dlg)
if ~exist(up.paths.annotations_folder, 'dir')
    mkdir(up.paths.annotations_folder)
end

%% waveform specifications
up.spikefs = 500;     % It is assumed in vortal_spike_conversion that all the spike waveforms are exported at this freq
up.reg = 1;
up.paworig = 2;
up.pawresample = 3;
up.resp = 4;
up.FS = [4, up.spikefs, 50, 62.5];

%% Operational specifications
up.DOSAVE = 1;                  % Save plots of filter characteristics?
up.DOPLOT = 0;
up.viewing_win_length = 32;     % duration of window to be viewed at once (in secs)
up.edge_prop = 0.3;             % proportion of window to show either side of the one currently being viewed.

%% Filtering specifications
up.filt_prop = 2;              % proportion of a second for which to apply the filter. (i.e. this * fs = order)
up.LP_CUTOFF = 60;                 % in rpm
up.time_offset = 1.1;          % The time (secs) by which the resp signal is delayed after the Paw signal.

%% Analysis specifications:
up.win_length = 32;            % window length for analysis in secs
up.true_BR_lims = 2;           % number of bpm either side at which ref BR estimates are considered reasonable.
up.mingap = 5;                 % If there is a gap in data of duration greater than this no. of secs, then this is considered to be a gap (anything less is ignored), for both numerics and waveforms

% duration of Tukey window taper in secs
up.paramSet.tukey_win_duration_taper = 2;

% Eliminate HFs (above resp freqs)
up.paramSet.elim_hf.Fpass = 1.2;  % in Hz
up.paramSet.elim_hf.Fstop = 0.8899;  % in Hz     (1.2 and 0.8899 provide a -3dB cutoff of 1 Hz)
up.paramSet.elim_hf.Dpass = 0.057501127785;
up.paramSet.elim_hf.Dstop = 0.01;

end

function [axis_h, h2] = annotation_plotter2(annotations_filepath, h, imp, grey_els, rel_els, s, subjects, up, win, NUMWINS, starttime)

%% ------------    annotation_plotter  ------------------------------
%
% Created by: Peter Charlton (17/10/2012)
%
% Purpose: To plot the Paw and impedance waveforms, and the annotated
%           breaths.
%
% Context: Some safety nets have been added to make it run ok even if the
% impedance signal isn't recorded.
%
% Inputs:
%
% Outputs:
%
% Files required:
%
% Called by:            Run_paw_verification
%
% ------------------------------------------------------------------------

ftsize = 16; lwidth = 2;
init_time = imp.t(rel_els(1));

load(annotations_filepath)
mean_imp = mean(imp.v(rel_els));
h3 = plot(imp.t(rel_els) - init_time, zeros(1,length(rel_els)), 'k', 'LineWidth', lwidth); hold on    % plot zero line on Paw plot
set(h3, 'HandleVisibility', 'off')
h1 = plot(imp.t(grey_els) - init_time, imp.v(grey_els) - mean_imp, 'Color', [0.5 0.5 0.5], 'LineWidth', lwidth); hold on,    % plot Paw (grey)
set(h1, 'HandleVisibility', 'off')
h2 = plot(imp.t(rel_els) - init_time, imp.v(rel_els) - mean_imp, 'LineWidth', lwidth);    % plot Paw
set(h2, 'HandleVisibility', 'off')
xlim([imp.t(grey_els(1)) - init_time, imp.t(grey_els(end)) - init_time])

%% Calculate and plot mix signal
if up.mode == 2
    % plot annotations which are relevant to this plot:
    %breathels = find(pk_anns.t>=imp.t(grey_els(1)) & pk_anns.t<=imp.t(grey_els(end)));
    %plot(pk_anns.t(breathels) - init_time,pk_anns.v(breathels)-mean_imp, 'ro','LineWidth',4)
    breathels = find(pk_anns.t>=imp.t(rel_els(1)) & pk_anns.t<=imp.t(rel_els(end)));
    plot(pk_anns.t(breathels) - init_time,pk_anns.v(breathels), 'ro','LineWidth',4)
    
    title(['Please annotate the peaks (red)  -  ' num2str(win) ' of ' num2str(NUMWINS) ' windows for subject ' num2str(s) ' of ' num2str(subjects(end))], 'FontSize', ftsize)
    ylabel('Impedance', 'FontSize', ftsize), xlabel('Time [s]', 'FontSize', ftsize);
    ylabel('Imepdance (above line = inhalation)', 'FontSize', ftsize)
    set(gca,'XTick',((imp.t(rel_els(1)) - init_time):5:ceil((imp.t(rel_els(end)) - init_time))), 'YTick', [])
    set(gca, 'FontSize', ftsize)
else
    title(['Please annotate high quality (left) or low quality (right)  -  ' num2str(win) ' of ' num2str(NUMWINS) ' windows for subject ' num2str(s) ' of ' num2str(subjects(end))], 'FontSize', ftsize)
    ylabel('Impedance', 'FontSize', ftsize), xlabel('Time [s]', 'FontSize', ftsize);
    ylabel('Imepdance (above line = inhalation)', 'FontSize', ftsize)
    set(gca,'XTick',((imp.t(rel_els(1)) - init_time):5:ceil((imp.t(rel_els(end)) - init_time))), 'YTick', [])
    set(gca, 'FontSize', ftsize)
    if find(qual.t == starttime)
        rel_qual = qual.v(qual.t == starttime);
        if rel_qual
            annotation('textbox',[0.75, 0.2, 0.1,0.1],'String','High Quality', 'Color', 'b', 'FontSize', ftsize*3, 'LineStyle', 'None')
        else
            annotation('textbox',[0.75, 0.2, 0.1,0.1],'String','Low Quality', 'Color', 'r', 'FontSize', ftsize*3, 'LineStyle', 'None')
        end
    end
    if exist('pk_anns', 'var')
        % plot annotations which are relevant to this plot:
        breathels = find(pk_anns.t>=imp.t(rel_els(1)) & pk_anns.t<=imp.t(rel_els(end)));
        plot(pk_anns.t(breathels) - init_time,pk_anns.v(breathels), 'ro','LineWidth',4)
    end
end

axis_h = gca;

end

function [pk_anns, TRS]= Click_CallBack2(h,e, a, axis_h, h2, up, annotations_filepath, start_time)

switch get(ancestor(a,'figure'),'SelectionType')
    
    case 'normal' %left click
        point = get(a,'CurrentPoint');
        load(annotations_filepath);
        [pk_anns.t,inds]=unique([pk_anns.t; point(1,1)+start_time]);
        pk_anns.v=[pk_anns.v; point(1,2)]; pk_anns.v = pk_anns.v(inds);
        [pk_anns.t, inds] = sort(pk_anns.t);
        pk_anns.v = pk_anns.v(inds);
        save(annotations_filepath, 'pk_anns');
        
    case 'alt'  % right click
        point = get(a,'CurrentPoint');
        load(annotations_filepath);
        [pk_anns.t,inds]=unique([pk_anns.t; point(1,1)]);
        pk_anns.v=[pk_anns.v; point(1,2)]; pk_anns.v = pk_anns.v(inds);
        [pk_anns.t, inds] = sort(pk_anns.t);
        pk_anns.v = pk_anns.v(inds);
        save(annotations_filepath, 'pk_anns');
        
    case 'extend' % right click whilst holding down shift
        load(annotations_filepath)
        pk_anns.t = pk_anns.t(1:(end-1));
        pk_anns.v = pk_anns.v(1:(end-1));
        save(annotations_filepath, 'pk_anns');
        
end

cla(axis_h)
plot(pk_anns.t-start_time,pk_anns.v, 'ro','LineWidth',4)

end

function paw_verification_script_checker(subjects, data, up)
% this function is used to check that the results can indeed be plotted
% back on the original waveform:
%% initialisation
fprintf('\n\n\n\n\n~~~~~~~  Reviewing Annotations  ~~~~~~~\n');

%% find out the file details
operator=input('Please type the first and last initial,\nof the operator whose annotations\nyou''d like to review, and press enter:\n','s');
question=input('\nWould you like to start from the \nfirst subject? (y or n, followed by enter)\n', 's');
if strcmp(question, 'n') == 1
    start_subject = input('\nWhich subject?\n', 's');
    start_subject = str2double(start_subject);
else
    start_subject = 1;
end

start_subject_el = find(subjects == start_subject);
if isempty(start_subject_el)
    fprintf('\nThis subject doesn''t exist. Try all over again.');
else
    subjects = subjects(start_subject_el:end);
    for s = subjects(1:end)
        
        %% Load Impedance
        imp.fs = data(s).ref.resp_sig.imp.fs;
        imp.v = data(s).ref.resp_sig.imp.v;
        imp.t = [0:(length(imp.v)-1)]/imp.fs;
        
        %% Process Impedance
        imp = lpf_to_exclude_resp(imp, up);
        imp.v = -1*imp.v;
        
        %% Load annotations
        annotations_filepath = [up.paths.annotations, 'RRLISTEN' num2str(s), '_' operator '-imp_an.mat'];
        up_copy = up;
        load(annotations_filepath);
        up = up_copy; clear up_copy
        
        %% plot waveforms and annotations
        duration_of_signal = imp.t(end) - imp.t(1);
        NUMWINS = floor(duration_of_signal / up.win_length); clear duration_of_signal
        edge_time = up.edge_prop*up.win_length;
        for win = 1:NUMWINS
            % Setup
            starttime = imp.t(1)+edge_time+(win-1)*up.win_length;
            endtime = starttime + up.win_length;
            rel_els = find(imp.t >= starttime & imp.t < endtime);
            grey_els = find(imp.t >= (starttime - edge_time) & imp.t < (endtime + edge_time));
            scrsz = get(0,'ScreenSize');
            
            % Plot
            h=figure('OuterPosition',[1 1 scrsz(3) scrsz(4)]);
            [axis_h, h2] = annotation_plotter2(annotations_filepath, h, imp, grey_els, rel_els, s, subjects, up, win, NUMWINS, starttime);
            
            % close plot after a pause
            pause
            close all
        end
    end
end

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
% norm_cutoff_freq = 0.3998;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(Fs/2);

% Prepare signal
s_dt=detrend(s);
s_filt = filtfilt(AMfilter.numerator, 1, s_dt);
end