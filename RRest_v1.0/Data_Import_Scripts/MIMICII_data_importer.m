function MIMICII_data_importer
% MIMICII_data_importer imports data from the MIMIC II database for use
% with RRest, a toolbox of respiratory rate algorithms.
%
%               MIMICII_data_importer
%
%	Inputs:
%       none
%                       just specify the relevant paths in the
%                       "universal_parameters" function below.
%
%	Outputs:
%       a single file containing all the data is written to the path
%       specified by "up.paths.analysis_path".
%
%   Requirements:
%       This requires the WFDB Toolbox, which can be downloaded from:
%           https://physionet.org/physiotools/matlab/wfdb-app-matlab/
%       data download is performed automatically using the script below.
%           
%   Further Information:
%       This version of the MIMICII_data_importer is provided to facilitate
%       reproduction of the analysis performed in:
%           Charlton P.H. et al., "Waveform Analysis to Estimate
%           Respiratory Rate" [In Press]
%       Further information on this study can be obtained at:
%           http://peterhcharlton.github.io/RRest/waveform_analysis.html
%       In addition, further information on RRest, including future
%       versions, can be obtained at:
%           http://peterhcharlton.github.io/RRest/index.html
%
%   Comments, Questions, Criticisms, Feedback, Contributions:
%       See: http://peterhcharlton.github.io/RRest/contributions.html
%
%   Version:
%       v.1 - published on 23rd Feb 2016 by Peter Charlton
%
%   Licence:
%       please see the accompanying file named "LICENSE"
%

%% Setup
up = universal_parameters;

%% Download data
download_data(up);

%% Identify patient stays
pt_stays = identify_pt_stays(up); 

%% Extract data from files on computer
data = extract_from_mimic_ii(pt_stays, up);

%% Save data in common format
save_data_in_common_format(data, up);

%% Identify sub-groups of neonates and adults
% This isn't necessary to obtain the data, but creates a plot used in the
% book chapter.
results = sub_group_identification(up);

end

function up = universal_parameters

fprintf('\n -- Setting up Universal Parameters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PARAMETERS TO BE SPECIFIED %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the root data directory (where the data will be stored)
up.paths.data_root = 'C:\Documents\Data\mimicii\';
up.paths.data_save_folder = 'C:\Documents\Data\';

% Specify the web address of the data to be downloaded
rel_database = 'mimic2wdb';
up.paths.database_dir = ['http://physionet.org/physiobank/database/', rel_database];

% Please note that your system may require the slashes in file paths to be
% of the opposite direction. In which case, change the following:
slash = '\';

% if you want .eps illustrations, then do as follows:
up.eps_figs = 0;  % set this to 1
if up.eps_figs
    % download 'export_fig' from:
    % http://uk.mathworks.com/matlabcentral/fileexchange/23629-export-fig
    % add the path of your download:
    export_fig_dir_path = 'C:\Users\pc13\Documents\GitHub\phd\Tools\Other Scripts\export_fig\';
    addpath(export_fig_dir_path)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   OTHER PARAMETERS   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   (don't need editing)   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extraction definitions
up.extraction.period = 600;    % period of waveform data to extract (in s)
up.extraction.rel_sigs = {'II', 'PLETH', 'RESP'};    % required signals
up.extraction.rel_nums = {'HR', 'PULSE', 'RESP'};    % required numerics
up.no_pt_stays = 100;

% database definitions
up.mimic_db.version = 3;
up.mimic_db.part = 0;

% paths
up.paths.file_location = [up.paths.data_root, 'physionet.org\physiobank\database\', rel_database, '\', num2str(up.mimic_db.version), num2str(up.mimic_db.part), '\'];
up.paths.analysis_path = [up.paths.data_root, 'Analysis_files', slash 'Data_for_Analysis', slash ];
up.paths.plots = [up.paths.data_root, 'Analysis_files', slash 'Results', slash 'Figures', slash ];
up.paths.records_filepath = [up.paths.analysis_path, 'RECORDS'];

% url paths
up.paths.url.records = [up.paths.database_dir, '/', num2str(up.mimic_db.version), num2str(up.mimic_db.part), '/', 'RECORDS'];

% Add all functions within the directory
addpath(genpath(fileparts(mfilename('fullpath'))));

% check that all the required folders exist. If not, create them:
req_folders = {up.paths.data_root, up.paths.analysis_path, up.paths.plots};
displayed_msg = 0;
for req_folder_no = 1 : length(req_folders)
    curr_folder = req_folders{req_folder_no};
    if ~exist(curr_folder, 'dir')
        if displayed_msg == 0
            fprintf('\n -- Creating folder in which to store data and analysis')
            displayed_msg = 1;
        end
        mkdir(curr_folder)
    end
end

% Check that the WFDB Matlab Toolbox has been installed
if ~exist('getWfdbClass', 'file')
    error('Couldn''t find the WFDB Matlab Toolbox. Please install as described at the top of the file.')
end

end

function download_data(up)

fprintf('\n -- Downloading Data')

%% set current directory to that of wget:
cd(up.paths.data_root)

%% Find out what records are available:
% download list of records
outfilename = websave(up.paths.records_filepath,up.paths.url.records);

% form list of records from downloaded file
records = load(up.paths.records_filepath);

%% Download each record
rel_records = records(1:up.no_pt_stays);

for rec_no = 1 : length(rel_records)
    
    % specify record
    record = rel_records(rec_no);
    
    % specify download location
    down_loc = [up.paths.database_dir, '/', num2str(up.mimic_db.version), num2str(up.mimic_db.part), '/', num2str(record), '/', 'RECORDS'];
    
    % set save location
    save_folder = [up.paths.file_location, num2str(record), '\'];
    save_loc = [save_folder, 'RECORDS'];
    if ~exist(save_folder, 'dir')
        mkdir(save_folder)
    end
    
    % download RECORDS file for this record
    outfilename = websave(save_loc,down_loc);
    
    % find out what individual files are in this record:
    fileID = fopen(save_loc);
    file_names = textscan(fileID, '%s %*[^\n]'); file_names = file_names{1};
    fclose(fileID);
    
    for file_no = 1 : length(file_names)
        
        % download header file
        curr_file_name = [file_names{file_no}, '.hea'];
        down_loc = [up.paths.database_dir, '/', num2str(up.mimic_db.version), num2str(up.mimic_db.part), '/', num2str(record), '/', curr_file_name];
        save_loc = [save_folder, curr_file_name];
        outfilename = websave(save_loc,down_loc);
        
        if ~strcmp(file_names{file_no}, num2str(record))
            % also download data file
            curr_file_name = [file_names{file_no}, '.dat'];
            down_loc = [up.paths.database_dir, '/', num2str(up.mimic_db.version), num2str(up.mimic_db.part), '/', num2str(record), '/', curr_file_name];
            save_loc = [save_folder, curr_file_name];
            outfilename = websave(save_loc,down_loc);
        end
    end
    
end

end

function pt_stays = identify_pt_stays(up)

fprintf('\n -- Identifying Patient Stays')

% Folders (each one is a patient stay)
dirs = dir(up.paths.file_location);

elim_els = ones(length(dirs),1);
for s = 1:length(dirs)
    if length(dirs(s).name) < 3
        elim_els(s) = 0;
    end
end

dirs = dirs(logical(elim_els'));

pt_stays = extractfield(dirs, 'name');


end

function data = extract_from_mimic_ii(pt_stays, up)
%
% Based on RDSAMP, part of the WFDB Toolbox. Further details are available
% at:
%       http://physionet.org/physiotools/matlab/wfdb-app-matlab/
%
% Reads MIMIC II data and outputs it as a structure of .t and .v values.

fprintf('\n -- Extracting Data from MIMIC II files')

if up.no_pt_stays > length(pt_stays)
    warning('There are fewer ICU stay datasets than requested')
end
no_pt_stays = min([up.no_pt_stays, length(pt_stays)]);

extraction_stats.no_pt_stays = no_pt_stays;
extraction_stats.no_records_searched = 0;
extraction_stats.no_records_without_signals = 0;
extraction_stats.no_records_without_num = 0;
extraction_stats.no_records_with_sigs_and_num = 0;

for pt_stay_no = 1 : no_pt_stays
    
    % set current dir
    dir_path = [up.paths.file_location, pt_stays{pt_stay_no}, '\'];
    cd(dir_path)
    
    % check to see if there are any files
    temp = dir([dir_path, pt_stays{pt_stay_no}, '_0*.hea']);
    if ~isempty(temp)
        records = extractfield(temp, 'name'); clear temp
    else
        extraction_stats.no_records_searched = extraction_stats.no_records_searched+1;
        extraction_stats.no_records_without_signals = extraction_stats.no_records_without_signals+1;
        continue
    end
    clear dir_path
    
    
    remaining_time_to_extract = up.extraction.period;
    
    for record_no = 1 : length(records)
        
        % setup java
        persistent javaWfdbExec config
        if(isempty(javaWfdbExec))
            [javaWfdbExec,config]=getWfdbClass('rdsamp');
        end
        
        %Remove file extension
        record_str = records{record_no}(1:(end-4));
        
        % Extract file information (including which signals are present)
        [siginfo,Fs]=extract_file_info(record_str);
        
        % Identify relevant signals
        sig_inds = 1:siginfo.no_sigs;
        sig_log = zeros(length(sig_inds),1);
        for sig_no = 1 : length(sig_inds)
            if sum(strcmp(up.extraction.rel_sigs, siginfo.sig_names(sig_no)))
                sig_log(sig_no) = 1;
            end
        end
        
        % skip if this doesn't have both the required signals
        if sum(sig_log) ~= length(up.extraction.rel_sigs)
            clear current_waves extracted_data Fs ListCapacity siginfo signalList javaWfdbExec config wfdb_argument
            continue
        end
        
        % create wfdb argument
        N0 = 1;
        wfdb_argument={'-r',record_str,'-Ps','-f',['s' num2str(N0-1)]};
        
        % Make command to only extract the required signals
        signalList = sig_inds(logical(sig_log));
        if(~isempty(signalList))
            wfdb_argument{end+1}='-s ';
            %-1 is necessary because WFDB is 0 based indexed.
            for sInd=1:length(signalList)
                wfdb_argument{end+1}=[num2str(signalList(sInd)-1)];
            end
        end
        
        % determine how many samples to extract from this file
        remaining_samples_to_extract = remaining_time_to_extract*Fs;
        N=siginfo.no_samps;
        if remaining_samples_to_extract < N
            N = remaining_samples_to_extract;
        end
        remaining_time_to_extract = remaining_time_to_extract - (N/Fs);
        
        % skip if we don't require any more data
        if N == 0
            clear current_waves extracted_data Fs ListCapacity siginfo signalList javaWfdbExec config wfdb_argument
            break
        end
        
        % additional bits of wfdb_argument
        wfdb_argument{end+1}='-t';
        wfdb_argument{end+1}=['s' num2str(N)];
        
        % extract data
        extracted_data=javaWfdbExec.execToDoubleArray(wfdb_argument);
        clear wfdb_argument
        
        %% Re-format data for output
        
        waves.fs = siginfo.fs;
        
        % store this new bit of the wave
        current_waves = fieldnames(waves);
        signalNames = siginfo.sig_names(signalList);
        for wave_no = 1 : length(up.extraction.rel_sigs)
            rel_sig_ind = find(strcmp(signalNames, up.extraction.rel_sigs{wave_no}));
            if ~strcmp(current_waves, up.extraction.rel_sigs{wave_no})
                eval(['waves.' up.extraction.rel_sigs{wave_no} '.t = extracted_data(:, 1);']);
                eval(['waves.' up.extraction.rel_sigs{wave_no} '.v = extracted_data(:, rel_sig_ind+1);']);
            else
                eval(['current_last_time = waves.' up.extraction.rel_sigs{wave_no} '.t(end);']);
                eval(['waves.' up.extraction.rel_sigs{wave_no} '.t = [waves.' up.extraction.rel_sigs{wave_no} '.t; current_last_time + (1/Fs) + extracted_data(:, 1)];']);
                eval(['waves.' up.extraction.rel_sigs{wave_no} '.v = [waves.' up.extraction.rel_sigs{wave_no} '.v; extracted_data(:, rel_sig_ind+1)];']);
            end
        end
        clear current_waves extracted_data Fs ListCapacity siginfo signalList javaWfdbExec config
        
    end
    
    
    %% Extract corresponding numerics
    
    if exist('waves', 'var')
        
        % setup java
        persistent javaWfdbExec_num config_num
        if(isempty(javaWfdbExec_num))
            [javaWfdbExec_num,config_num]=getWfdbClass('rdsamp');
        end
        
        filepath = [pt_stays{pt_stay_no}, 'n'];
        
        if exist([up.paths.file_location, pt_stays{pt_stay_no}, '\', filepath, '.hea'], 'file')
            
            % Extract file information (including which signals are present)
            [numinfo,Fs]=extract_file_info(filepath);
            loc = numinfo.loc;
            
            % Identify relevant numerics
            num_inds = 1:numinfo.no_sigs;
            num_log = zeros(length(num_inds),1);
            for num_no = 1 : length(num_inds)
                if sum(strcmp(up.extraction.rel_nums, numinfo.sig_names(num_no)))
                    num_log(num_no) = 1;
                end
            end
            
            % skip if this doesn't have both the required signals
            if sum(num_log) == length(up.extraction.rel_nums)
                
                numerics.fs = numinfo.fs;
                
                % create wfdb argument
                N0 = 1;
                N = up.extraction.period*numerics.fs;
                wfdb_argument={'-r',filepath,'-Ps','-f',['s' num2str(N0-1)]};
                
                % Make command to only extract the required signals
                numList = num_inds(logical(num_log));
                if(~isempty(numList))
                    wfdb_argument{end+1}='-s ';
                    %-1 is necessary because WFDB is 0 based indexed.
                    for sInd=1:length(numList)
                        wfdb_argument{end+1}=[num2str(numList(sInd)-1)];
                    end
                end
                
                % Additional wfdb_argument bits
                wfdb_argument{end+1}='-t';
                wfdb_argument{end+1}=['s' num2str(N)];
                
                % extract data
                extracted_data=javaWfdbExec_num.execToDoubleArray(wfdb_argument);
                clear wfdb_argument
                
                %% Re-format data for output
                
                % store this new bit of the wave
                current_numerics = fieldnames(numerics);
                numNames = numinfo.sig_names(numList);
                for wave_no = 1 : length(up.extraction.rel_nums)
                    rel_num_ind = find(strcmp(numNames, up.extraction.rel_nums{wave_no}));
                    if ~strcmp(current_numerics, up.extraction.rel_nums{wave_no})
                        eval(['numerics.' up.extraction.rel_nums{wave_no} '.t = extracted_data(:, 1);']);
                        eval(['numerics.' up.extraction.rel_nums{wave_no} '.v = extracted_data(:, rel_num_ind+1);']);
                    else
                        eval(['current_last_time = numerics.' up.extraction.rel_nums{wave_no} '.t(end);']);
                        eval(['numerics.' up.extraction.rel_nums{wave_no} '.t = [numerics.' up.extraction.rel_nums{wave_no} '.t; current_last_time + (1/Fs) + extracted_data(:, 1)];']);
                        eval(['numerics.' up.extraction.rel_nums{wave_no} '.v = [numerics.' up.extraction.rel_nums{wave_no} '.v; extracted_data(:, rel_num_ind+1)];']);
                    end
                end
                
            end
            
        end
        
        clear current_numerics extracted_data Fs ListCapacity numinfo numList javaWfdbExec_num config_num
        
    end
    
    
    %% Find extraction statistics
    extraction_stats.no_records_searched = extraction_stats.no_records_searched+1;
    if ~exist('numerics', 'var')
        extraction_stats.no_records_without_num = extraction_stats.no_records_without_num+1;
    end
    if ~exist('waves', 'var')
        extraction_stats.no_records_without_signals = extraction_stats.no_records_without_signals+1;
    end
    
    %% Store data for this patient stay in the output
    if exist('waves', 'var') && exist('numerics', 'var')
        extraction_stats.no_records_with_sigs_and_num = extraction_stats.no_records_with_sigs_and_num+1;
        data{extraction_stats.no_records_with_sigs_and_num}.loc = loc;
        data{extraction_stats.no_records_with_sigs_and_num}.ID = pt_stays{pt_stay_no};
        data{extraction_stats.no_records_with_sigs_and_num}.waves = waves;
        data{extraction_stats.no_records_with_sigs_and_num}.numerics = numerics;
        clear waves numerics
    end
    
end

%% Find out how many records had the required duration signals
total_records = length(data);
recs_w_req_duration = zeros(total_records,1);
for s = 1  : length(data)
    if sum(strcmp(fieldnames(data{1,s}), 'waves'))
        if length(data{1,s}.waves.II.t) == up.extraction.period*data{1,s}.waves.fs ...
            && length(data{1,s}.waves.PLETH.t) == up.extraction.period*data{1,s}.waves.fs ...
            && length(data{1,s}.waves.RESP.t) == up.extraction.period*data{1,s}.waves.fs
            recs_w_req_duration(s) = 1;            
        end
    end
end

%% Only output those records with the required duration
extraction_stats.no_records_with_req_duration = sum(recs_w_req_duration);
data = data(logical(recs_w_req_duration));

%% Create Pie Chart of extraction stats

% bar chart
ftsize = 14;
figure('Position', [200, 200, 800, 450])
labels = {'Signals or numerics absent', 'Shorter than 10 mins', 'Included in analysis'};
piedata(3) = extraction_stats.no_records_with_req_duration;
piedata(2) = extraction_stats.no_records_with_sigs_and_num - extraction_stats.no_records_with_req_duration;
piedata(1) = extraction_stats.no_records_searched - extraction_stats.no_records_with_sigs_and_num;
p = pie(piedata);
for s = 2:2: length(p)
    p(s).FontSize = ftsize;
end
set(gca, 'FontSize', ftsize)
legend(labels, 'Location', 'southeastoutside')
title('Reasons for exclusion of records', 'FontSize', ftsize)
% save
save_name = 'record_exclusion';
savepath = [up.paths.plots, save_name, '.png'];
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[1 1 8 4.5]);
print(gcf, savepath, '-dpng')
close

end

function [headers, Fs] = extract_file_info(filepath)

%% Obtain info on signals

fileID = fopen([filepath, '.hea']);
header_lines = textscan(fileID,'%s','Delimiter','\n'); header_lines = header_lines{1,1};
fclose all;

% Extract the relevant header information
temp.top_line = textscan(header_lines{1}, '%s', 'delimiter', sprintf(' ')); temp.top_line = temp.top_line{1,1};  % adapted from http://uk.mathworks.com/matlabcentral/newsreader/view_thread/249016
headers.no_sigs = str2double(temp.top_line{2});
headers.fs = str2double(temp.top_line{3});
headers.no_samps = str2double(temp.top_line{4});

% extract signal names
headers.sig_names = cell(0);
for sig_no = 1 : headers.no_sigs
    temp.curr_line = textscan(header_lines{1+sig_no}, '%s', 'delimiter', sprintf(' ')); temp.curr_line = temp.curr_line{1,1};  % adapted from http://uk.mathworks.com/matlabcentral/newsreader/view_thread/249016
    headers.sig_names = [headers.sig_names; temp.curr_line{end}];    
end

% extract location
loc_line = header_lines{end};
headers.loc = loc_line(13:end);

% putting in a format for rdsamp

Fs = headers.fs;

end

function save_data_in_common_format(old_data, up)

fprintf('\n -- Saving data in appropriate format')

for subj_el = 1:length(old_data)
    
    SID = old_data{1,subj_el}.ID;
    loc = old_data{1,subj_el}.loc;
    
    % insert group name
    if strcmp(loc, 'nicu')
        data(1,subj_el).group = 'neonate';
    else
        data(1,subj_el).group = 'adult';
    end
    
    % Insert fixed params
    
    data(1,subj_el).fix.id = SID;
    data(1,subj_el).fix.loc = loc;
    data(1,subj_el).fix.ventilation = 'unknown';
    data(1,subj_el).fix.recording_conditions = 'critical care';
    
    % insert PPG signal
    ppg = old_data{1,subj_el}.waves.PLETH;
    data(1,subj_el).ppg.v = ppg.v;
    clear ppg
    data(1,subj_el).ppg.fs = old_data{1,subj_el}.waves.fs;
    % insert EKG signal
    ekg = old_data{1,subj_el}.waves.II;
    data(1,subj_el).ekg.v = ekg.v;
    clear ekg
    data(1,subj_el).ekg.fs = old_data{1,subj_el}.waves.fs;
    % insert RESP signal
    resp = old_data{1,subj_el}.waves.RESP;
    data(1,subj_el).ref.resp_sig.imp.v = resp.v;
    clear resp
    data(1,subj_el).ref.resp_sig.imp.fs = old_data{1,subj_el}.waves.fs;
    data(1,subj_el).ref.resp_sig.imp.method = 'thoracic impedance measured using clinical monitor';
    
    % insert numerics
    data(1,subj_el).ref.params.rr.t = old_data{1,subj_el}.numerics.RESP.t;
    data(1,subj_el).ref.params.rr.v = old_data{1,subj_el}.numerics.RESP.v;
    data(1,subj_el).ref.params.rr.method = 'derived from thoracic impedance measured using clinical monitor';
    data(1,subj_el).ref.params.rr.units.t = 's';
    data(1,subj_el).ref.params.rr.units.v = 'breaths/min';
    data(1,subj_el).ref.params.hr.t = old_data{1,subj_el}.numerics.HR.t;
    data(1,subj_el).ref.params.hr.v = old_data{1,subj_el}.numerics.HR.v;
    data(1,subj_el).ref.params.hr.method = 'derived from ecg measured using clinical monitor';
    data(1,subj_el).ref.params.hr.units.t = 's';
    data(1,subj_el).ref.params.hr.units.v = 'beats/min';
    data(1,subj_el).ref.params.pr.t = old_data{1,subj_el}.numerics.PULSE.t;
    data(1,subj_el).ref.params.pr.v = old_data{1,subj_el}.numerics.PULSE.v;
    data(1,subj_el).ref.params.pr.method = 'derived from ppg measured using clinical monitor';
    data(1,subj_el).ref.params.pr.units.t = 's';
    data(1,subj_el).ref.params.pr.units.v = 'beats/min';
    
end

% Save to file
save([up.paths.data_save_folder, 'mimicii_data'], 'data')

end

function results = sub_group_identification(up)

%% Load data
load([up.paths.data_save_folder, 'mimicii_data'], 'data')

%% Split into sub-groups and identify data for each
sub_group = zeros(length(data),1);   % zero for neonates, one for adults
n_rrs = []; a_rrs = [];
for s = 1 : length(data)
    if strcmp(data(1,s).group, 'neonate')
        n_rrs = [n_rrs; data(1,s).ref.params.rr.v];
    else
        a_rrs = [a_rrs; data(1,s).ref.params.rr.v];
        sub_group(s) = 1;
    end
end
results.sub_group = sub_group;

%% Plot histograms

% generate data of prop of rrs in each category:
xbins1 = -1:120;
rrs.a.tot = sum(~isnan(a_rrs));
rrs.a.bins = nan(length(xbins1),1);
rrs.n.tot = sum(~isnan(n_rrs));
rrs.n.bins = nan(length(xbins1),1);
for s = 1:length(xbins1)
    rel_rr_range = [(xbins1(s)-0.5),(xbins1(s)+0.5)];
    rel_a_rrs = a_rrs(a_rrs>rel_rr_range(1) & a_rrs<=rel_rr_range(2));
    rel_n_rrs = n_rrs(n_rrs>rel_rr_range(1) & n_rrs<=rel_rr_range(2));
    rrs.a.bins(s) = length(rel_a_rrs)/rrs.a.tot;
    rrs.n.bins(s) = length(rel_n_rrs)/rrs.n.tot;
end

% bar chart
ftsize = 14;
paper_size = [9, 3.5];
figure('Position', [200, 200, 100*paper_size(1), 100*paper_size(2)])
b = bar(xbins1, [rrs.a.bins, rrs.n.bins], 1.2);
b(1).FaceColor = [0,0,1];
b(1).EdgeColor = [0,0,1];
b(2).FaceColor = [1, 0.6, 0.6];
b(2).EdgeColor = [1, 0.6, 0.6];
xlim([-1 xbins1(end)])
ylim([0 1.1*max([rrs.a.bins; rrs.n.bins])])
xlabel('Impedance RR numeric [bpm]', 'FontSize', ftsize)
ylabel('Proportion of Measurements', 'FontSize', ftsize)
set(gca, 'YTick', [], 'FontSize', ftsize)
legend(b, {'Adults', 'Neonates'})

%% save histograms
save_name = 'adults_vs_neonates';
savepath = [up.paths.plots, save_name];
PrintFigs(gcf, paper_size, savepath, up)

end

function PrintFigs(h, paper_size, savepath, up)
set(h,'PaperUnits','inches');
set(h,'PaperSize', [paper_size(1), paper_size(2)]);
set(h,'PaperPosition',[0 0 paper_size(1) paper_size(2)]);
print(h,'-dpdf',savepath)
print(h,'-dpng',savepath)

if up.eps_figs
    export_fig(savepath, '-eps')
end
close all
end