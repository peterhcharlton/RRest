function ImP_SQI_mimic_data_importer
% ImP_SQI_mimic_data_importer imports data from the MIMIC II database with
% which to assess the performance of the impedance pneumography SQI, and
% outputs the dataset in Matlab format.
%
%               ImP_SQI_mimic_data_importer
%
%	Inputs:
%       Just specify the relevant paths in the "universal_parameters"
%       function below. No other inputs are required as this script
%       downloads the required data files from PhysioNet automatically. 
%
%	Outputs:
%       A Matlab file containing the collated data, called "mimic_imp_sqi_data.mat"
%
%   The WFDB Toolbox:
%       This requires the WFDB Toolbox, which can be downloaded from:
%           https://archive.physionet.org/physiotools/matlab/wfdb-app-matlab/
%       I installed it using the 'Quick Start' instructions at:
%           https://archive.physionet.org/physiotools/matlab/wfdb-app-matlab/#quick-start
%               (including modifying the download path to be https://archive.physionet.org ... )
%       The Toolbox is now available at:
%           https://physionet.org/content/wfdb-swig-matlab/1.0.0/
%           
%   Further Information:
%       This version of the ImP_SQI_mimic_data_importer is provided to facilitate
%       reproduction of the MIMIC subset used in:
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
up = universal_parameters;

%% Create ReadMe and License files
create_readme_and_license_files(up);

%% Download data
download_data(up);

%% Extract data from files on computer
extract_from_mimic_ii(up);

%% Add in manual annotations
add_manual_annotations(up);

% %% Save data in RRest format
% if up.output.rrest_format, save_data_in_rrest_format(up); end
% 
% %% Save data in CSV format
% if up.output.csv_format, save_data_in_csv_format(up); end
% 
% %% Save data in WFDB format
% if up.output.csv_format, save_data_in_wfdb_format(up); end

fprintf(['\n\n -- Finished importing the ImP SQI dataset. Data saved at:\n\n       ', up.paths.data_root, '\n\n'])

end

function up = universal_parameters

fprintf('\n -- Setting up Universal Parameters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PARAMETERS TO BE SPECIFIED %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the root data directory (where the data will be stored)
up.paths.data_root = '/Users/petercharlton/Documents/Data/Imp_SQI_mimic_dataset_test/';
fprintf('\n   - The data will be stored at %s', up.paths.data_root);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   OTHER PARAMETERS   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   (don't need editing)   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

provide_details_of_script;

% Specify the web address of the data to be downloaded
rel_database = 'mimic3wdb/matched';
rel_database_part = 'p00';
up.paths.database_dir_root = ['https://archive.physionet.org/physiobank/database/', rel_database, '/'];
up.paths.database_dir = [up.paths.database_dir_root, rel_database_part];
up.paths.analysis_path = [up.paths.data_root, 'Analysis_files', filesep, 'Data_for_Analysis', filesep ];
up.paths.records_filepath = [up.paths.analysis_path, 'RECORDS'];
up.paths.download_details_file = [up.paths.analysis_path, 'download_details'];
up.no_pt_stays = 100;

% extraction definitions
up.extraction.rel_sigs = {'II', 'PLETH', 'RESP'};    % required signals
up.extraction.rel_nums = {'HR', 'PULSE', 'RESP', 'SpO2'};    % required numerics
up.fs = 125;  % Hz
up.req_data_time = 0*60;  % data taken from 0 minutes into record (secs)
up.req_data_duration = 60*60;  % segments of data extracted (secs)
up.dataset_name = 'imp_sqi_mimic';
up.csv.new_line = '\n';
up.units.signals = {'II', 'AVR', 'V',  'III', 'I',  'MCL1', 'MCL', 'PLETH', 'RESP', 'ABP',  'CVP',  'ART',   'P1',  'UAP',  'PAP'};
up.units.units   = {'mV', 'mV',  'mV', 'mV',  'mV', 'mV',   'mV',  'NU',    'pm'  , 'mmHg', 'mmHg', 'mmHg', 'mmHg', 'mmHg', 'mmHg'};
up.units.numerics    = {'HR',  'PULSE',  'SpO2', 'RESP', 'ABPSys', 'ABPDias', 'ABPMean', 'NBPSys', 'MBPDia', 'MBPMean'};
up.units.num_units   = {'bpm', 'bpm',    '%',    'pm',   'mmHg',   'mmHg',    'mmHg',    'mmHg',   'mmHg',   'mmHg'   };

% paths
up.paths.file_location = [up.paths.data_root, 'physionet.org', filesep, 'physiobank' ,filesep 'database', filesep];
up.paths.extracted_data_file = [up.paths.data_root, 'mimic_imp_sqi_data.mat'];

% url paths
up.paths.url.records = [up.paths.database_dir_root, 'RECORDS'];

% Add all functions within the directory
addpath(genpath(fileparts(mfilename('fullpath'))));

% check that all the required folders exist. If not, create them:
if ~exist(up.paths.data_root, 'dir')
    fprintf('\n The specified data folder does not exist \n - please either: \n 1) select a different location (speifying it in "universal_parameters"), or\n 2) create a folder at the specified location: %s', up.paths.data_root);
    error('Specified folder does not exist')
end
req_folders = {up.paths.analysis_path};
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

% Type(s) of output (currently not used)
up.output.wfdb_format = true; % Save in PhysioNet's waveform database format
up.output.csv_format = true;  % Save as comma-separated-value format
up.output.rrest_format = true;% Save in a format suitable for use with the RRest Toolbox of algorithms

% Check that the WFDB Matlab Toolbox has been installed
if ~exist('getWfdbClass', 'file')
    error('Couldn''t find the WFDB Matlab Toolbox. Please install as described at the top of the file.')
end

end

function provide_details_of_script

licence_details = ['\n\n ImP_SQI_data_importer', ...
    '\n This program comes with ABSOLUTELY NO WARRANTY', ... 
    '\n and is available under the GNU public license.', ...
    '\n For details see the accompanying LICENSE.txt file.\n\n'];

fprintf(licence_details)

end

function download_data(up)

fprintf('\n -- Downloading Data:')

%% set current directory to that of wget:
cd(up.paths.data_root)

%% Find out what records are available:
fprintf(' List of records, ');
% download list of records
if ~exist(up.paths.records_filepath, 'file')
    outfilename = websave(up.paths.records_filepath,up.paths.url.records);
end

% form list of patient stays from downloaded file
pt_stays = importdata(up.paths.records_filepath);

%% Download each record
file_counter = 0;

for pt_stay_no = 1 : length(pt_stays)
    
    % specify patient stay
    pt_stay = pt_stays{pt_stay_no};
    temp = strfind(pt_stay, '/');
    database_subset = pt_stay(1:temp(1)-1);
    pt_no = pt_stay(temp(1)+1:end-1); clear temp
    
    % set save location
    save_folder = [up.paths.file_location, pt_stay];
    if ~exist(save_folder, 'dir')
        mkdir(save_folder)
    end
    
    % download records for this pt stay
    pt_stay_records_url = [up.paths.database_dir_root, pt_stay, 'RECORDS'];
    save_path = [save_folder,'RECORDS'];
    if ~exist(save_path, 'file')
        outfilename = websave(save_path, pt_stay_records_url);
    end
    pt_stay_records = importdata(save_path);
    
    % download overall header file
    pt_stay_segment = pt_stay_records{1};
    curr_file_name = [pt_stay_segment, '.hea'];
    overall_header_file = [up.paths.database_dir_root, pt_stay, curr_file_name];
    overall_header_file_save_loc = [save_folder, curr_file_name];
    if ~exist(overall_header_file_save_loc, 'file')
        outfilename = websave(overall_header_file_save_loc,overall_header_file);
    end
    
    % find out what individual files are in this pt stay:
    fileID = fopen(overall_header_file_save_loc);
    individual_files.name = {};
    individual_files.duration = []; individual_files.cum_duration_before = []; individual_files.cum_duration_after = [];
    line_no = 0;
    while ~feof(fileID)
        curr_text_line = fgetl(fileID); line_no = line_no+1;
        if line_no == 1
            temp = strfind(curr_text_line, ':');
            individual_files.start_time = curr_text_line(temp(1)-2:temp(1)+9);
        elseif ~isempty(strfind(curr_text_line, '_')) & length(strfind(curr_text_line, ' '))==1 & isempty(strfind(curr_text_line, 'layout'))
            individual_files.name{end+1} = curr_text_line(1:strfind(curr_text_line, '_')+4);
            individual_files.duration(end+1) = str2double(curr_text_line(strfind(curr_text_line, '_')+6:end));
        end
        clear curr_text_line
    end
    clear line_no
    fclose(fileID); clear fileID
    
    % Find out which files contain sufficient data
    rel_file_els = find(individual_files.duration >= (up.req_data_duration*up.fs+1));
    rel_file_names = individual_files.name(rel_file_els);
    
    % cycle through each file containing sufficient data
    for file_no = 1 : length(rel_file_names)
        curr_rel_file = rel_file_names{file_no};
        
        % download data header file
        curr_file_name = [curr_rel_file, '.hea'];
        down_loc = [up.paths.database_dir_root, pt_stay, curr_file_name];
        save_loc = [save_folder, curr_file_name];
        if ~exist(save_loc, 'file')
            outfilename = websave(save_loc,down_loc);
        end
        data_header_file = save_loc;
        
        % identify the signals in this file
        fileID = fopen(save_loc);
        line_no = 0;
        signals = [];
        while ~feof(fileID)
            curr_text_line = fgetl(fileID); line_no = line_no+1;
            if line_no == 1
                temp = strfind(curr_text_line, ':');
                curr_file.start_time = curr_text_line(temp(1)-2:end);
            elseif contains(curr_text_line, '.dat')
                temp = strfind(curr_text_line, ' ');
                signals{end+1} = curr_text_line(temp(end)+1:end);
            end
            clear curr_text_line
        end
        clear line_no
        fclose(fileID); clear fileID
        
        % check to see whether this file contains the required signals
        if sum(strcmp(signals, 'II')) && sum(strcmp(signals, 'PLETH')) && sum(strcmp(signals, 'RESP'))
            
            % store details
            file_counter = file_counter+1;
            fprintf([' ' num2str(file_counter)]);
            download.source.waves_header_file{file_counter} = down_loc;
            download.source.overall_header_file{file_counter} = overall_header_file;
            download.saved.pt_stay{file_counter} = pt_stay;
            download.saved.pt_no{file_counter} = pt_no;
            download.saved.pt_stay_segment{file_counter} = pt_stay_segment;
            download.saved.waves_header_file{file_counter} = data_header_file;            
            download.saved.root_folder{file_counter} = save_folder;
            download.saved.overall_header_file{file_counter} = overall_header_file_save_loc;
            
            % download data file
            curr_file_name = [curr_rel_file, '.dat'];
            download.saved.file_name{file_counter} = curr_file_name;
            down_loc = [up.paths.database_dir_root, pt_stay, curr_file_name];
            download.source.waves_data_file{file_counter} = down_loc;
            download.source.waves_data_filename{file_counter} = curr_file_name;
            save_loc = [save_folder, curr_file_name];
            if ~exist(save_loc, 'file')
                outfilename = websave(save_loc,down_loc);
            end
            download.saved.waves_data_file{file_counter} = save_loc;
            clear outfilename save_loc down_loc curr_file_name
            
            % download numerics header file
            curr_file_name = [pt_stay_segment, 'n.hea'];
            down_loc = [up.paths.database_dir_root, pt_stay, curr_file_name];
            save_loc = [save_folder, curr_file_name];
            if ~exist(save_loc, 'file')
                outfilename = websave(save_loc,down_loc);
                download.source.numerics_header_file{file_counter} = down_loc;
            end
            download.saved.numerics_header_file{file_counter} = save_loc;
            clear outfilename save_loc down_loc curr_file_name
            
            % download numerics data file
            temp = strfind(curr_rel_file, '_');
            curr_file_name = [curr_rel_file(1:temp-1), 'n.dat']; clear temp
            down_loc = [up.paths.database_dir_root, pt_stay, curr_file_name];
            download.source.numerics_data_file{file_counter} = down_loc;
            save_loc = [save_folder, curr_file_name];
            if ~exist(save_loc, 'file')
                outfilename = websave(save_loc,down_loc);
            end
            download.saved.numerics_data_file{file_counter} = save_loc;
            clear outfilename save_loc down_loc curr_file_name
            
            % skip any remaining files for this subject
            break
            
        end
        
        clear outfilename save_loc down_loc curr_file_name
    end
    
    % finish if reached the desired number of pt stays with data
    if file_counter >= up.no_pt_stays
        break
    end
    
end

% Save the download details
save(up.paths.download_details_file, 'download');

end

function extract_from_mimic_ii(up)
%
% Based on RDSAMP, part of the WFDB Toolbox. Further details are available
% at:
%       http://physionet.org/physiotools/matlab/wfdb-app-matlab/
%
% Reads MIMIC II data and outputs it as a structure of .t and .v values.

fprintf('\n -- Extracting Data from MIMIC files: ')

% Load download details
load(up.paths.download_details_file);

for pt_stay_no = 1 : length(download.saved.pt_stay)
    fprintf([' ' num2str(pt_stay_no)]);
    
    % set current dir
    dir_path = download.saved.root_folder{pt_stay_no};
    cd(dir_path)
    
    % setup java
    % persistent javaWfdbExec config
    [javaWfdbExec, config] = deal([]);
    if(isempty(javaWfdbExec))
        [javaWfdbExec,config]=getWfdbClass('rdsamp');
    end
    
    % Remove file extension
    record_str = download.saved.file_name{pt_stay_no}(1:end-4);
    
    % Extract file information (including which signals are present)
    siginfo = extract_file_info(record_str);
    
    % create wfdb argument
    N0 = (siginfo.fs * up.req_data_time);
    wfdb_argument={'-r',record_str,'-Ps','-f',['s' num2str(N0)]}; clear record_str
    
    % calculate start time of numerics data
    [numerics_start, waveforms_start] = calc_start_times(download, pt_stay_no);
    
    % Make command to only extract the required signals
    wfdb_argument{end+1}='-s ';
    %-1 is necessary because WFDB is 0 based indexed.
    for sInd=1:siginfo.no_sigs
        wfdb_argument{end+1} = num2str(sInd-1);
    end
    clear sInd
    
    % determine how many samples to extract from this file
    N= N0 +(siginfo.fs * up.req_data_duration);
    
    % additional bits of wfdb_argument
    wfdb_argument{end+1}='-t';
    wfdb_argument{end+1}=['s' num2str(N+1)];
    
    % extract waveform data
    extracted_data=javaWfdbExec.execToDoubleArray(wfdb_argument);
    clear wfdb_argument
    
    % Store waveform data
    waves.fs = siginfo.fs;
    for wave_no = 1 : length(siginfo.sig_names)
        eval(['waves.' siginfo.sig_names{wave_no} '.t = extracted_data(:, 1)-(N0/125);']);
        eval(['waves.' siginfo.sig_names{wave_no} '.v = extracted_data(:, wave_no+1);']);
    end
    clear wave_no N N0 current_waves extracted_data Fs ListCapacity signalList javaWfdbExec config
    
    
    %% Extract corresponding numerics
        
    % Extract file information (including which signals are present)
    header_filepath = download.saved.numerics_header_file{pt_stay_no};
    do_numerics = 1;
    if ~isempty(header_filepath) && do_numerics
        header_filepath = header_filepath(1:end-4);
        numinfo = extract_file_info(header_filepath);
        
        % setup java
        persistent javaWfdbExec_num config_num
        if(isempty(javaWfdbExec_num))
            [javaWfdbExec_num,config_num]=getWfdbClass('rdsamp');
        end
        
        % Identify relevant numerics
        num_inds = 1:numinfo.no_sigs;
        num_log = zeros(length(num_inds),1);
        for num_no = 1 : length(num_inds)
            if sum(strcmp(up.extraction.rel_nums, numinfo.sig_names(num_no)))
                num_log(num_no) = 1;
            end
        end
        
        numerics.fs = numinfo.fs;
        
        % create wfdb argument
        if length(numinfo.start_time) == 5, numinfo.start_time = ['00:', numinfo.start_time]; end
%         start_time_diff = (datenum(['00/01/0000 ' numinfo.start_time], 'dd/mm/yyyy HH:MM:SS') - datenum(['00/01/0000 ' siginfo.start_time], 'dd/mm/yyyy HH:MM:SS')) *24*60*60;  % in secs
        start_time_diff = round((waveforms_start - numerics_start) *24*60*60);  % in secs
        N0 = (numerics.fs * up.req_data_time) + (start_time_diff*numerics.fs);
        clear waveforms_start numerics_start start_time_diff
        N = N0 +(numerics.fs * up.req_data_duration);
        [~, header_filepath, ~] = fileparts(header_filepath);
        wfdb_argument={'-r',header_filepath,'-Ps','-f',['s' num2str(N0)]}; clear header_filepath
        
        % Make command to only extract the required signals
        numList = num_inds(logical(num_log)); clear num_log num_inds
        if(~isempty(numList))
            wfdb_argument{end+1}='-s ';
            %-1 is necessary because WFDB is 0 based indexed.
            for sInd=1:length(numList)
                wfdb_argument{end+1}=[num2str(numList(sInd)-1)];
            end
            clear sInd
        end
        
        % Additional wfdb_argument bits
        wfdb_argument{end+1}='-t';
        wfdb_argument{end+1}=['s' num2str(N+1)];
        
        % extract numerics
        extracted_data = javaWfdbExec_num.execToDoubleArray(wfdb_argument);
        clear wfdb_argument
        
        % store the numerics
        numNames = numinfo.sig_names(numList);
        for num_no = 1 : length(up.extraction.rel_nums)
            rel_num_ind = find(strcmp(numNames, up.extraction.rel_nums{num_no}));
            eval(['numerics.' up.extraction.rel_nums{num_no} '.t = round(extracted_data(:, 1)-N0);']);
            eval(['numerics.' up.extraction.rel_nums{num_no} '.v = extracted_data(:, rel_num_ind+1);']);
            clear rel_num_ind
        end
        data{pt_stay_no}.numerics = numerics;
        
        clear numerics num_no numinfo siginfo current_numerics extracted_data numNames ListCapacity numinfo numList javaWfdbExec_num config_num N N0
        
    end
    
    %% Extract fixed data
    
    fileID = fopen(download.saved.overall_header_file{pt_stay_no});
    while ~feof(fileID)
        curr_text_line = fgetl(fileID);
        if ~isempty(strfind(curr_text_line, 'Location'))
            fix.loc = curr_text_line(strfind(curr_text_line, ': ')+2:end);
        end
        if ~isempty(strfind(curr_text_line, 'age'))
            fix.age = curr_text_line(strfind(curr_text_line, 'age>: ')+6:strfind(curr_text_line, '<sex>: ')-2);
            if isempty(strfind(fix.age, '?')) && isempty(strfind(fix.age, '+'))
                fix.age = str2double(fix.age);
            elseif strcmp(fix.age, '??')
                fix.age = nan;
            end
        end
        if ~isempty(strfind(curr_text_line, 'sex'))
            fix.sex = curr_text_line(strfind(curr_text_line, 'sex>: ')+6:end);
        end
        clear curr_text_line
    end
    clear line_no
    fclose(fileID); clear fileID
    
    %% Store data for this patient stay in the output
    data{pt_stay_no}.loc = fix.loc;
    if sum(strcmp(fieldnames(fix), 'age'))
        data{pt_stay_no}.age = fix.age;
    end
    if sum(strcmp(fieldnames(fix), 'sex'))
        data{pt_stay_no}.sex = fix.sex;
    end
    data{pt_stay_no}.ID = download.saved.pt_stay{pt_stay_no};
    data{pt_stay_no}.data_segment = download.saved.pt_stay_segment{pt_stay_no};
    data{pt_stay_no}.waves = waves;
    clear waves numerics fix
    
end

%% Convert to required variable for analysis

old_data = data; clear data

for s = 1:length(old_data)
    data(s).ref.resp_sig.imp.v = old_data{s}.waves.RESP.v;
    data(s).ref.resp_sig.imp.fs = old_data{s}.waves.fs;
    data(s).ref.params.hr.t = old_data{s}.numerics.HR.t;
    data(s).ref.params.hr.v = old_data{s}.numerics.HR.v;
    data(s).ref.params.hr.units.t = 'secs since beginning of signal';
    data(s).ref.params.hr.units.v = 'beats per minute';
    data(s).ref.params.pr.t = old_data{s}.numerics.PULSE.t;
    data(s).ref.params.pr.v = old_data{s}.numerics.PULSE.v;
    data(s).ref.params.pr.units.t = 'secs since beginning of signal';
    data(s).ref.params.pr.units.v = 'beats per minute';
    data(s).ref.params.rr.t = old_data{s}.numerics.RESP.t;
    data(s).ref.params.rr.v = old_data{s}.numerics.RESP.v;
    data(s).ref.params.rr.units.t = 'secs since beginning of signal';
    data(s).ref.params.rr.units.v = 'breaths per minute';
    data(s).ref.params.spo2.t = old_data{s}.numerics.SpO2.t;
    data(s).ref.params.spo2.v = old_data{s}.numerics.SpO2.v;
    data(s).ref.params.spo2.units.t = 'secs since beginning of signal';
    data(s).ref.params.spo2.units.v = '%';
    data(s).ekg.v = old_data{s}.waves.II.v;
    data(s).ekg.fs = old_data{s}.waves.fs;
    data(s).ppg.v = old_data{s}.waves.PLETH.v;
    data(s).ppg.fs = old_data{s}.waves.fs;
    data(s).fix.id = old_data{s}.ID;
    data(s).fix.loc = old_data{s}.loc;
    data(s).fix.source.waves = download.source.waves_data_file{s};
    data(s).fix.source.numerics = download.source.numerics_data_file{s};
    
end

%% Save extracted data
save(up.paths.extracted_data_file, 'data')

end

function [numerics_start, waveforms_start] = calc_start_times(download, pt_stay_no)

% Identify relevant header files
all_waves_header_file = download.saved.overall_header_file{pt_stay_no};
specific_waves_header_file = download.saved.waves_header_file{pt_stay_no};
numerics_header_file = download.saved.numerics_header_file{pt_stay_no};

% Find start time of all waveform data
fid = fopen(all_waves_header_file);
lines = textscan(fid,'%s','Delimiter','\n'); header_line = lines{1,1}{1};
fclose(fid);
colons = strfind(header_line, ':');
timings.all_waves.start_time = header_line(colons(1)-2:colons(2)+6);
slashes = strfind(header_line, '/'); slashes = slashes(slashes>colons(2));
timings.all_waves.start_date = header_line(slashes(1)-2:slashes(2)+4);
timings.all_waves.deb = datenum([timings.all_waves.start_date, ' ', timings.all_waves.start_time], 'dd/mm/yyyy HH:MM:SS.FFF');
clear slashes colons fid lines

% Find start time of the waveform data in this file
fid = fopen(all_waves_header_file);
lines = textscan(fid,'%s','Delimiter','\n'); lines = lines{1,1};
fclose(fid);
% First line is a header line:
lines = lines(2:end);
% Exclude next line if it's a layout line
if contains(lines{1}, 'layout 0')
    lines = lines(2:end);
end
% Exclude last line if it's a location line
if strcmp(lines{end}(1), '#')
    lines = lines(1:end-1);
end
% Identify duration (and cumulative duration) of each waveform file
for line_no = 1 : length(lines)
    space_loc = strfind(lines{line_no}, ' ');
    line_data_duration(line_no,1) = str2double(lines{line_no}(space_loc+1:end));
    line_data_name{line_no,1} = lines{line_no}(1:space_loc-1);
end
clear fid lines line_no space_loc
line_cum_duration = [0;cumsum(line_data_duration)]; % in ms
% Find start time of the waveform file of interest
rel_line = find(strcmp(line_data_name, download.source.waves_data_filename{pt_stay_no}(1:end-4)));
timings.rel_waves.deb = timings.all_waves.deb + (line_cum_duration(rel_line)/(125*60*60*24));

% % Find start time of waveforms data
% fid = fopen(specific_waves_header_file);
% lines = textscan(fid,'%s','Delimiter','\n'); header_line = lines{1,1}{1};
% fclose(fid);
% colons = strfind(header_line, ':');
% timings.rel_waves.start_time = header_line(colons(1)-2:colons(2)+6);
% slashes = strfind(header_line, '/'); slashes = slashes(slashes>colons(2));
% timings.rel_waves.start_date = header_line(slashes(1)-2:slashes(2)+4);
% timings.rel_waves.deb = datenum([timings.rel_waves.start_date, ' ', timings.rel_waves.start_time], 'dd/mm/yyyy HH:MM:SS');
% clear slashes colons fid lines

% Find start time of all numerics data
fid = fopen(numerics_header_file);
lines = textscan(fid,'%s','Delimiter','\n'); header_line = lines{1,1}{1};
fclose(fid);
colons = strfind(header_line, ':');
timings.all_numerics.start_time = header_line(colons(1)-2:colons(2)+6);
slashes = strfind(header_line, '/'); slashes = slashes(slashes>colons(2));
timings.all_numerics.start_date = header_line(slashes(1)-2:slashes(2)+4);
timings.all_numerics.deb = datenum([timings.all_numerics.start_date, ' ', timings.all_numerics.start_time], 'dd/mm/yyyy HH:MM:SS');
clear slashes colons fid lines

numerics_start = timings.all_numerics.deb;
waveforms_start = timings.rel_waves.deb;

end

function headers = extract_file_info(filepath)

%% Obtain info on signals

fileID = fopen([filepath, '.hea']);
header_lines = textscan(fileID,'%s','Delimiter','\n'); header_lines = header_lines{1,1};
fclose all;

% Extract the relevant header information
temp.top_line = textscan(header_lines{1}, '%s', 'delimiter', sprintf(' ')); temp.top_line = temp.top_line{1,1};  % adapted from http://uk.mathworks.com/matlabcentral/newsreader/view_thread/249016
headers.no_sigs = str2double(temp.top_line{2});
headers.fs = str2double(temp.top_line{3});
headers.no_samps = str2double(temp.top_line{4});
rel_line = find(~cellfun(@isempty, strfind(temp.top_line, ':')));
headers.start_time = temp.top_line{rel_line}; 
clear rel_line

% extract signal names
headers.sig_names = cell(0);
for sig_no = 1 : headers.no_sigs
    temp.curr_line = textscan(header_lines{1+sig_no}, '%s', 'delimiter', sprintf(' ')); temp.curr_line = temp.curr_line{1,1};  % adapted from http://uk.mathworks.com/matlabcentral/newsreader/view_thread/249016
    headers.sig_names = [headers.sig_names; temp.curr_line{end}];    
end

% extract location
loc_line = header_lines{end};
headers.loc = loc_line(13:end);

end

function save_data_in_rrest_format(up)

fprintf('\n -- Saving data in RRest toolbox format')

temp = dir([up.paths.individual_files_folder, '*.mat']);
created_files = extractfield(temp, 'name'); clear temp

for subj_el = 1:length(created_files)
    
    % Load this subject's data
    subj_data = load([up.paths.individual_files_folder, created_files{subj_el}]);
    
    % Insert fixed params
    data(1,subj_el).fix.id = created_files{subj_el}(10:15);
    data(1,subj_el).fix.data_segment = created_files{subj_el}(17:end);
    data(1,subj_el).fix.loc = subj_data.location;
    data(1,subj_el).fix.ventilation = 'unknown';
    data(1,subj_el).fix.recording_conditions = 'critical care';
    data(1,subj_el).fix.source = subj_data.download_source;
    
    % insert PPG signal
    rel_col = strcmp(subj_data.header, 'PLETH');
    data(1,subj_el).ppg.v = subj_data.data(:,rel_col);
    data(1,subj_el).ppg.fs = subj_data.fs(rel_col);
    data(1,subj_el).ppg.method = 'measured using clinical monitor';
    % insert EKG signal
    rel_col = strcmp(subj_data.header, 'II');
    data(1,subj_el).ekg.v = subj_data.data(:,rel_col);
    data(1,subj_el).ekg.fs = subj_data.fs(rel_col);
    data(1,subj_el).ekg.method = 'measured using clinical monitor';
    % insert RESP signal
    rel_col = strcmp(subj_data.header, 'RESP');
    data(1,subj_el).ref.resp_sig.imp.v = subj_data.data(:,rel_col);
    data(1,subj_el).ref.resp_sig.imp.fs = subj_data.fs(rel_col);
    data(1,subj_el).ref.resp_sig.imp.method = 'thoracic impedance measured using clinical monitor';
    % insert annotations
    data(1,subj_el).ref.breaths.ann1 = subj_data.ann.resp.breaths.x1;
    data(1,subj_el).ref.breaths.ann2 = subj_data.ann.resp.breaths.x2;
    data(1,subj_el).ref.breaths.method = 'two independent sets of breath annotations using the impedance respiratory signal';
    
    % insert numerics
    if sum(strcmp(fieldnames(subj_data), 'numerics'))
        data(1,subj_el).ref.params.rr.t = subj_data.numerics.RESP.t;
        data(1,subj_el).ref.params.rr.v = subj_data.numerics.RESP.v;
        data(1,subj_el).ref.params.rr.method = 'derived from thoracic impedance measured using clinical monitor';
        data(1,subj_el).ref.params.rr.units.t = 's';
        data(1,subj_el).ref.params.rr.units.v = 'breaths/min';
        data(1,subj_el).ref.params.hr.t = subj_data.numerics.HR.t;
        data(1,subj_el).ref.params.hr.v = subj_data.numerics.HR.v;
        data(1,subj_el).ref.params.hr.method = 'derived from ecg measured using clinical monitor';
        data(1,subj_el).ref.params.hr.units.t = 's';
        data(1,subj_el).ref.params.hr.units.v = 'beats/min';
        data(1,subj_el).ref.params.pr.t = subj_data.numerics.PULSE.t;
        data(1,subj_el).ref.params.pr.v = subj_data.numerics.PULSE.v;
        data(1,subj_el).ref.params.pr.method = 'derived from ppg measured using clinical monitor';
        data(1,subj_el).ref.params.pr.units.t = 's';
        data(1,subj_el).ref.params.pr.units.v = 'beats/min';
        data(1,subj_el).ref.params.spo2.t = subj_data.numerics.SpO2.t;
        data(1,subj_el).ref.params.spo2.v = subj_data.numerics.SpO2.v;
        data(1,subj_el).ref.params.spo2.method = 'measured using clinical monitor';
        data(1,subj_el).ref.params.spo2.units.t = 's';
        data(1,subj_el).ref.params.spo2.units.v = '%';
    end
    
end

% Save to file
savepath = [up.paths.data_root, up.dataset_name, '_data'];
save(savepath, 'data')

fprintf(['\n --- Dataset saved to: ' savepath])

end

function save_data_in_csv_format(up)

fprintf('\n -- Converting to CSV format');

% lower case of dataset name to keep to usual format of PhysioNet files:
dataset_name = up.dataset_name;

% Create folder to save data in
data_type = 'csv';
save_folder = [up.paths.data_root, dataset_name, '_', data_type, filesep];
if ~exist(save_folder, 'dir')
    mkdir(save_folder)
end
cd(save_folder)

% Load data
load(up.paths.extracted_data_file);

% cycle through each subject
no_subjs = length(data);
file_names = cell(no_subjs,1);
for subj_no = 1 : no_subjs
    
    % setup command
    file_names{subj_no} = [up.dataset_name, '_', num2str(subj_no, '%02.f')];
    fs = data{subj_no}.waves.fs;
    sig_names = fieldnames(data{subj_no}.waves); sig_names = sig_names(~strcmp(sig_names, 'fs'));
    num_names = fieldnames(data{subj_no}.numerics); num_names = num_names(~strcmp(num_names, 'fs'));
    
    % create a file of fixed variables for this subject
    curr_filename = [file_names{subj_no}, '_Fix.txt'];
    fid = fopen(curr_filename, 'wt');
    
    signals_line = 'Time [s],';
    data_mat = data{subj_no}.waves.II.t;
    for sig_no = 1 : length(sig_names)
        signals_line = [signals_line, ' ', sig_names{sig_no}];
        eval(['curr_sig = data{subj_no}.waves.' sig_names{sig_no} ';'])
        data_mat = [data_mat, curr_sig.v];
        % check whether time vector is consistent
        if ~isequal(curr_sig.t, data{subj_no}.waves.II.t)
            error('Check this')
        end
        
        if sig_no < length(sig_names)
            signals_line = [signals_line, ';'];
        end
    end
    clear sig_no curr_sig
    
    numerics_line = 'Time [s],';
    numerics_mat = data{subj_no}.numerics.HR.t;
    numerics_fs = data{subj_no}.numerics.fs;
    for num_no = 1 : length(num_names)
        numerics_line = [numerics_line, ' ', num_names{num_no}];
        eval(['curr_num = data{subj_no}.numerics.' num_names{num_no} ';'])
        numerics_mat = [numerics_mat, curr_num.v];
        % check whether time vector is consistent
        if ~isequal(curr_num.t, data{subj_no}.numerics.HR.t)
            error('Check this')
        end
        
        if num_no < length(num_names)
            numerics_line = [numerics_line, ';'];
        end
    end
    clear num_no curr_num
    
    % create a text file containing the fixed variables for this subject
    text_to_write = [file_names{subj_no}, up.csv.new_line, ...
        'Signals:' signals_line(10:end), up.csv.new_line, ...
        'Signals sampling frequency: ' num2str(fs) ' Hz', up.csv.new_line, ...
        'Numerics:' numerics_line(10:end), up.csv.new_line, ...
        'Numerics sampling frequency: ' num2str(numerics_fs) ' Hz', up.csv.new_line, ...
        'Age: ' num2str(data{subj_no}.age), up.csv.new_line, ...
        'Gender: ' data{subj_no}.sex, up.csv.new_line, ...
        'Location: ' data{subj_no}.loc, up.csv.new_line, ...
        'MIMIC II matched wdb ID: ' data{subj_no}.ID, up.csv.new_line, ...
        'MIMIC II data segment: ' data{subj_no}.data_segment, up.csv.new_line, ...
        'Original overall header file: ' data{subj_no}.download_source.overall_header_file, up.csv.new_line, ...
        'Original signals header file: ' data{subj_no}.download_source.waves_header_file, up.csv.new_line, ...
        'Original signals data file: ' data{subj_no}.download_source.waves_data_file, up.csv.new_line, ...
        'Original numerics header file: ' data{subj_no}.download_source.numerics_header_file, up.csv.new_line, ...
        'Original numerics data file: ' data{subj_no}.download_source.numerics_data_file, up.csv.new_line, ...
        ];
    fprintf(fid, text_to_write); 
    fclose(fid);
    clear text_to_write fid
    
    % create a csv file of signals for this subject
    curr_filename = [file_names{subj_no}, '_Signals.csv'];
    fid = fopen(curr_filename,'w');
    signals_line = strrep(signals_line, ';', ',');
    fprintf(fid,[signals_line, up.csv.new_line]); clear signals_line
    fclose(fid); clear fid
    dlmwrite(curr_filename, data_mat, '-append'); clear data_mat curr_filename
    
    % create a csv file of numerics for this subject
    curr_filename = [file_names{subj_no}, '_Numerics.csv'];
    fid = fopen(curr_filename,'w');
    numerics_line = strrep(numerics_line, ';', ',');
    fprintf(fid,[numerics_line, up.csv.new_line]); clear numerics_line
    fclose(fid); clear fid
    dlmwrite(curr_filename, numerics_mat, '-append'); clear numerics_mat curr_filename
    
    % create a csv file of annotations for this subject
    curr_filename = [file_names{subj_no}, '_Breaths.csv'];
    fid = fopen(curr_filename,'w');
    header_line = 'breaths ann1 [signal sample no], breaths ann2 [signal sample no]';
    fprintf(fid,[header_line, up.csv.new_line]); clear header_line
    fclose(fid); clear fid
    length_diff = length(data{subj_no}.ann.resp.breaths.x1(:)) - length(data{subj_no}.ann.resp.breaths.x2(:));
    if length_diff > 0
        data{subj_no}.ann.resp.breaths.x2 = [data{subj_no}.ann.resp.breaths.x2, nan(1, length_diff)];
    elseif length_diff < 0
        data{subj_no}.ann.resp.breaths.x1 = [data{subj_no}.ann.resp.breaths.x1, nan(1, abs(length_diff))];
    end
    ann_mat = [data{subj_no}.ann.resp.breaths.x1(:), data{subj_no}.ann.resp.breaths.x2(:)];
    dlmwrite(curr_filename, ann_mat, '-append'); clear ann_mat curr_filename
    
    clear fs sig_names signals_line curr_filename numerics_line
    
end
clear subj_no

end

function save_data_in_wfdb_format(up)

fprintf('\n -- Converting to WFDB format');

% Create folder to save data in
data_type = 'wfdb';
save_folder = [up.paths.data_root, up.dataset_name, '_', data_type, filesep];
if ~exist(save_folder, 'dir')
    mkdir(save_folder)
end
cd(save_folder)

% Load data
load(up.paths.extracted_data_file);

% cycle through each subject
no_subjs = length(data);
file_names = cell(no_subjs,1);
for subj_no = 1 : no_subjs
    
    % setup command
    file_names{subj_no} = [up.dataset_name, sprintf('%.2d', subj_no)];
    fs = data{subj_no}.waves.fs;
    sig_names = fieldnames(data{subj_no}.waves); sig_names = sig_names(~strcmp(sig_names, 'fs'));
    num_names = fieldnames(data{subj_no}.numerics); num_names = num_names(~strcmp(num_names, 'fs'));
    
    
    % Extract signals data
    data_mat = [];
    for sig_no = 1 : length(sig_names)
        eval(['curr_sig = data{subj_no}.waves.' sig_names{sig_no} ';'])
        data_mat = [data_mat, curr_sig.v];
        % check whether time vector is consistent
        if ~isequal(curr_sig.t, data{subj_no}.waves.II.t)
            error('Check this')
        end
    end
    clear sig_no curr_sig
    
    % Create signals file
    descrip = ['<age>: ' num2str(data{subj_no}.age) ' <sex>: ' data{subj_no}.sex ' <location>: ' data{subj_no}.loc ' <source>: https://physionet.org/physiobank/database/mimic2wdb/matched/' data{subj_no}.ID '/ <data segment>: ' data{subj_no}.data_segment '<modifications>: this is not a verbatim copy of an original file. Please see the accompanying LICENSE.txt file for further details.'];
    units = '';
    for sig_no = 1 : length(sig_names)
        signal_names(sig_no,1) = {[sig_names{sig_no} ', ' ]};
        rel_unit = up.units.units{strcmp(up.units.signals, sig_names{sig_no})};
        units = [units, rel_unit, '/']; clear rel_unit
    end
    clear sig_no
    units = units(1:end-1);
    % convert this subject's data into WFDB format
    mat2wfdb(data_mat, file_names{subj_no}, fs, [], units, descrip, [], signal_names);
    clear data_mat fs units descrip
    
    % Extract numerics data
    numerics_mat = [];
    for num_no = 1 : length(num_names)
        eval(['curr_num = data{subj_no}.numerics.' num_names{num_no} ';'])
        numerics_mat = [numerics_mat, curr_num.v];
        % check whether time vector is consistent
        if ~isequal(curr_num.t, data{subj_no}.numerics.HR.t)
            error('Check this')
        end
    end
    clear num_no curr_num
    
    % Create numerics file
    descrip = ['<age>: ' num2str(data{subj_no}.age) ' <sex>: ' data{subj_no}.sex ' <location>: ' data{subj_no}.loc ' <source>: ' data{subj_no}.download_source.numerics_data_file ' <modifications>: this is not a verbatim copy of an original file. Please see the accompanying LICENSE.txt file for further details.'];
    units = '';
    for num_no = 1 : length(num_names)
        numeric_names(num_no,1) = {[num_names{num_no} ', ' ]};
        rel_unit = up.units.num_units{strcmp(up.units.numerics, num_names{num_no})};
        units = [units, rel_unit, '/']; clear rel_unit
    end
    clear sig_no
    units = units(1:end-1);
    % convert this subject's data into WFDB format
    fs = data{subj_no}.numerics.fs;
    mat2wfdb(numerics_mat, [file_names{subj_no}, 'n'], fs, [], units, descrip, [], numeric_names);
    clear numerics_mat fs units descrip
    
    % create annotation files
    ann1 = data{subj_no}.ann.resp.breaths.x1;
    ann2 = data{subj_no}.ann.resp.breaths.x2;
    ann = [ann1(:); ann2(:)]; % col vector of samples of anns
    chan = find(~cellfun(@isempty, strfind(sig_names, 'RESP')))*ones(size(ann))-1;
    num = zeros(size(ann));
    comments = [repmat({'ann1'}, [length(ann1),1]); repmat({'ann2'}, [length(ann2),1])]; % ann1 or ann2 as char
    type = repmat('"', size(ann));
    subtype = zeros(size(num));
    wrann(file_names{subj_no},'breath',ann,type,subtype,chan,num,comments);
    
    clear sig_names ann1 ann2 ann chan num type subtype 
    
end
clear no_subjs subj_no

% create a RECORDS file
file_name = 'RECORDS';
file_path = [save_folder, file_name];
fid = fopen(file_path, 'w');
formatSpec = '%s\r\n';
for line_no = 1 : length(file_names)
    fprintf(fid, formatSpec, file_names{line_no});
end

% create an ANNOTATORS file
file_name = 'ANNOTATORS';
file_path = [save_folder, file_name];
fid = fopen(file_path, 'w');
fprintf(fid, 'ann1\tbreath annotations by annotator 1\n');
fprintf(fid, 'ann2\tbreath annotations by annotator 2');

% create a DBS file
file_name = 'DBS';
file_path = [save_folder, file_name];
fid = fopen(file_path, 'w');
fprintf(fid, [up.dataset_name, '\t', 'BIDMC ECG, PPG and Respiration Dataset']);

fclose all;

end


function create_readme_and_license_files(up)

fprintf('\n -- Saving ReadMe and License files')

%% License file
licence_file_path = [up.paths.data_root, 'LICENSE'];
fid = fopen(licence_file_path, 'wt');
fprintf(fid, ['The mimic_imp_sqi_data file produced by this script contains modified versions of PhysioNet files. The modified files are distributed under the terms for distribution of modified data files listed at: https://archive.physionet.org/copying.shtml#data .',...
    '\n\n The following modifications have been made:', ...
    '\n\n - Waveform, numerics and fixed data (such as age) were downloaded from the PhysioNet MIMIC III matched waveform database. 1-hour segments of data were extracted for each of the 100 ICU stays included in the dataset.',...
    '\n\n - The data were imported into Matlab. In addition, manually annotated breaths are provided for some 32 second windows of most subjects'' impedance pneumography signals, and estimated respiratory rates (calculated from these manual annotations) are provided.', ...
    '\n\n Further details of the dataset are provided in the accompanying ReadMe file.']);
fclose all;

readme_file_path = [up.paths.data_root, 'README'];
fid = fopen(readme_file_path, 'wt');

%% ReadMe file

readme_text = [...
'Impedance Pneumography Signal Quality Index Dataset',...
'\n=============',...
'\n\nTable of Contents',...
'\n=================',...
'\n\n  * [The Dataset](#the-dataset)',...
'\n  * [Revisions](#revisions)',...
'\n  * [Description](#description)',...
'\n\n# The Dataset',...
'\n\nThe Impedance Pneumography Signal Quality Index Dataset was first reported in the following publication:',...
'\n\nCharlton, P.H. et al. An impedance pneumography signal quality index: design, assessment and application to respiratory rate monitoring, [under review], 2020',...
'\n\nIn this publication it was used to assess the performance of a signal quality index for the impedance pneumography signal. The dataset was extracted from the much larger [MIMIC III matched waveform Database](https://archive.physionet.org/physiobank/database/mimic3wdb/matched/), which was acquired from critically-ill patients during hospital care at the Beth Israel Deaconess Medical Centre (Boston, MA, USA). The 100 recordings within the dataset, each of 1-hour duration, each contain:',...
'\n  * Physiological signals: the impedance pneumography respiratory signal (ImP), photoplethysmogram (PPG) and electrocardiogram (ECG). These are sampled at 125 Hz.',...
'\n  * Physiological parameters, such as the heart rate (HR), pulse rate (PR), respiratory rate (RR), and blood oxygen saturation level (SpO2). These are sampled at 1 Hz.',...
'\n  * Fixed parameters, such as age, sex, ward and patient ID',...
'\n  * Manual annotations of breaths, and RRs calculated from these: A single annotator manually annotated individual breaths in approximately five 32 sec segments of each recording using the impedance respiratory signal.',...
'\n\nThe dataset is distributed in Matlab (r) format.',...
'\n\nFor more information about the dataset, please contact author at:',...
'\npc657@medschl.cam.ac.uk .',...
'\n\n# Revisions',...
'\n\nR1: 	2017-Aug-06 		initial release',...
'\n\n## Description',...
'\n### Matlab (r) Format',...
'\nThe *mimic_imp_sqi_data.mat* file contains the following subset of the dataset in a single Matlab (r) variable named *data*. The following are provided for each of the 100 recordings:',...
'\n* *ekg*:   Lead II ECG signal. Each signal is provided in a structure, where the *v* field denotes the signal values, and *fs* is the sampling frequency.',...
'\n* *ppg*:   Photoplethysmogram signal',...
'\n* *ref.resp_sig.imp*:  Impedance respiratory signal',...
'\n* *ref.ann*:  Manual annotations of breaths provided by a single annotator. A vector of sample numbers is provided (*breath_els*), which correspond to the signal sample numbers. Breaths were only annotated in certain segments of the signal, rather than the whole signal: the start and end samples of each segment in which breaths were annotated are provided in *win_start_els* and *win_end_els*. RRs calculated from the manual annotations are provided (*rr*, one per window).',...
'\n* *ref.params*:  Physiological parameters: *rr* (respiratory rate, derived by the monitor from the impedance signal, breaths per minute), *hr* (heart rate, derived from the ECG, beats per minute), *pr* (pulse rate, derived from the PPG, beats per minute), *spo2* (blood oxygen saturation level, %%).',...
'\n* *fix*: A structure of fixed variables: *id* (the MIMIC III matched waveform database subject ID and recording identifier), *loc* (the ward location), and download source (specified for both the waves and numerics).',...
'\n\nFurther details on the contents of each file are provided [here](https://archive.physionet.org/tutorials/creating-records.shtml).'];

fprintf(fid, readme_text);
fclose all;

end

function add_manual_annotations(up)

add_manual_annotations = 1;

if ~add_manual_annotations, return, else clear add_manual_annotations, end

fprintf('\n -- Adding in Manual Annotations')

% load extracted data file
load(up.paths.extracted_data_file)

% specify folder containing manual annotations files
ann_folder = '/Users/petercharlton/Documents/Data/Imp_SQI_mimic_dataset_test/2019_annotations/';
manual_annotator = 'PC';
temp = dir([ann_folder, '*', manual_annotator, '_an.mat']);
ann_files = extractfield(temp, 'name'); clear temp

do_check = 0;
% import annotations for each manually annotated subject
for subj_no = 1 : length(data)
    
    % create fields in which to store annotations
    data(subj_no).ref.ann.win_start_els = [];
    data(subj_no).ref.ann.win_end_els = [];
    data(subj_no).ref.ann.breath_els = [];
    data(subj_no).ref.ann.rr = [];
    
    % identify candidate annotations file for this subject
    annotations_filename = ['ImP_SQI' num2str(subj_no) '_' manual_annotator '_an.mat'];
    annotations_filepath = [ann_folder, annotations_filename];
    
    % skip if no annotations for this subject
    if ~exist(annotations_filepath)
        continue
    end
    
    % load annotations (which are stored in 'pk_anns' and 'qual')
    temp = load(annotations_filepath);
    
    % identify windows
    imp.v = data(subj_no).ref.resp_sig.imp.v; imp.fs = 125;
    imp.t = [0:(length(imp.v)-1)]/imp.fs;
    up.paramSet.winLeng = 32;     % the duration of each window in secs (no overlap)
    win_starts = imp.t(1):up.paramSet.winLeng:imp.t(end);
    win_ends = win_starts+up.paramSet.winLeng;
    
    for win_no = 1 : length(win_starts)
        
        % Identify annotations within this window
        if sum(strcmp(fieldnames(temp), 'pk_anns'))
            win_start = win_starts(win_no);
            win_end = win_ends(win_no);
            rel_ann_els = temp.pk_anns.t >= win_start & temp.pk_anns.t <= win_end;
            if sum(rel_ann_els) > 1
                first_breath = min(temp.pk_anns.t(rel_ann_els));
                last_breath = max(temp.pk_anns.t(rel_ann_els));
                no_breaths = sum(rel_ann_els)-1;
                rr_ann = 60/((last_breath-first_breath)/no_breaths);
                
                data(subj_no).ref.ann.win_start_els = [data(subj_no).ref.ann.win_start_els; win_start*125+1];
                data(subj_no).ref.ann.win_end_els = [data(subj_no).ref.ann.win_end_els; win_end*125+1];
                data(subj_no).ref.ann.breath_els = [data(subj_no).ref.ann.breath_els; round(temp.pk_anns.t(rel_ann_els)*125+1)];
                data(subj_no).ref.ann.rr = [data(subj_no).ref.ann.rr; rr_ann];
                
                % to check:
                if do_check
                    sig_els = (win_starts(win_no)*125+1):(win_ends(win_no)*125+1); plot(imp.t(sig_els), imp.v(sig_els)), hold on, rel_breath_els = data(subj_no).ref.ann.breath_els(data(subj_no).ref.ann.breath_els >= sig_els(1) & data(subj_no).ref.ann.breath_els <= sig_els(end)); plot(imp.t(rel_breath_els), imp.v(rel_breath_els), '*r')
                    close all
                end
                
                clear rr_ann no_breaths first_breath last_breath
            end
            clear win_start win_end rel_ann_els
        end
        
    end
    % to check
    if do_check
        plot(imp.t,imp.v), hold on, plot(imp.t(data(subj_no).ref.ann.breath_els), imp.v(data(subj_no).ref.ann.breath_els), '*r')
        close all
    end

    clear win_no temp imp win_starts win_ends
    
end

% Save new data file with annotations added
save(up.paths.extracted_data_file, 'data')

end