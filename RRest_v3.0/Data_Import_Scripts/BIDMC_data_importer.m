function BIDMC_data_importer
% BIDMC_data_importer imports data from the MIMIC II database, and
% constructs the BIDMC dataset in a range of formats.
%
%               BIDMC_data_importer
%
%	Inputs:
%       Just specify the relevant paths in the "universal_parameters"
%       function below. No other inputs are required as this script
%       downloads the required data files from PhysioNet automatically. 
%
%	Outputs:
%       The following dataset files:
%           -   WFDB (WaveForm DataBase) Format files
%           -   a Matlab Format file
%           -   CSV (Comma-Separated Value) files
%
%   Requirements:
%       This requires the WFDB Toolbox, which can be downloaded from:
%           https://physionet.org/physiotools/matlab/wfdb-app-matlab/
%           
%   Further Information:
%       This version of the BIDMC_data_importer is provided to facilitate
%       reproduction of the dataset used in:
%           Pimental M. A. F. et al., "Towards a Robust Estimation of
%           Respiratory Rate from Pulse Oximeters", IEEE TBME, vol. 64(8),
%           pp.1914-1923, 2016. DOI: http://doi.org/10.1109/TBME.2016.2613124
%       Further information on this study can be obtained at:
%           http://peterhcharlton.github.io/RRest/bidmc_dataset.html
%       In addition, further information on RRest, a toolbox of respiratory
%       algorithms which can be used with this dataset, can be obtained at:
%           http://peterhcharlton.github.io/RRest/index.html
%
%   Comments, Questions, Criticisms, Feedback, Contributions:
%       See: http://peterhcharlton.github.io/RRest/contributions.html
%
%   Version:
%       v.1 - published on 20th June 2018 by Peter H Charlton
%
%   Source:
%       This script has been adapted from 'MIMICII_data_importer.m', which
%       is one of the scripts in the RRest toolbox, available at:
%           https://github.com/peterhcharlton/RRest
%
%   Licence:
%       please see the accompanying file named "LICENSE"
%

%% Setup
up = universal_parameters;

%% Create ReadMe and License files
create_readme_and_license_files(up);

%% Download data
download_data(up);

%% Extract data from files on computer
extract_from_mimic_ii(up);

%% Save data in individual files
save_data_to_individual_files(up);

%% Check whether individual files are correct
% used to check this data extraction script
do_check = 0; if do_check, check_individual_files(up); end

%% Save data in RRest format
if up.output.rrest_format, save_data_in_rrest_format(up); end

%% Save data in CSV format
if up.output.csv_format, save_data_in_csv_format(up); end

%% Save data in WFDB format
if up.output.csv_format, save_data_in_wfdb_format(up); end

fprintf(['\n\n -- Finished importing the BIDMC dataset. Data saved at:\n\n       ', up.paths.data_root, '\n\n'])

end

function up = universal_parameters

fprintf('\n -- Setting up Universal Parameters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PARAMETERS TO BE SPECIFIED %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the root data directory (where the data will be stored)
up.paths.data_root = '/Users/petercharlton/Documents/Data/BIDMC_dataset/';

% Type(s) of output
up.output.wfdb_format = true; % Save in PhysioNet's waveform database format
up.output.csv_format = true;  % Save as comma-separated-value format
up.output.rrest_format = true;% Save in a format suitable for use with the RRest Toolbox of algorithms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   OTHER PARAMETERS   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   (don't need editing)   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

provide_details_of_script;

% Specify the web address of the data to be downloaded
rel_database = 'mimic2wdb';
rel_database_part = 'matched';
up.paths.database_dir = ['http://physionet.org/physiobank/database/', rel_database, '/', rel_database_part];

% Direction of slashes in file paths:
slash = filesep;

% extraction definitions
up.extraction.rel_sigs = {'II', 'PLETH', 'RESP'};    % required signals
up.extraction.rel_nums = {'HR', 'PULSE', 'RESP', 'SpO2'};    % required numerics
up.fs = 125;  % Hz
up.pt_id =           {'s01182',          's01241',          's01795',          's01840',          's01995',          's03386',          's03386',          's03386',          's03386',          's03516',          's04641','s06914','s07432','s07784','s08452','s08915','s08936','s09483','s09968','s11342','s11342','s11342','s11342','s12531','s12589','s13314','s13927','s15646','s17236','s17497','s17735','s17756','s19981','s20471','s20486','s22348','s24455','s25323','s25323','s26845','s26868','s29093','s29125','s29503','s29622','s29866','s30026','s30243','s31106','s31141','s31400','s32195','s32628'};
up.pt_data_segment = {'2688-03-25-23-14','3195-12-18-05-38','2503-07-03-17-43','3454-10-24-18-46','3123-12-10-15-38','2577-07-30-13-04','2577-08-02-01-41','2577-08-03-12-35','2577-08-07-20-40','3087-07-13-01-12','3112-08-01-14-02','2530-02-22-13-21','3359-07-25-15-04','3115-08-30-15-19','2862-12-24-16-10','3391-06-02-00-27','3043-07-05-18-01','3089-08-26-21-25','2564-01-26-21-54','2777-07-19-19-33','2777-09-25-23-35','2777-09-28-13-38','2777-09-30-12-41','3144-03-08-13-21','3365-11-16-23-47','2641-10-27-08-25','3107-03-28-07-43','3164-11-16-15-38','2899-06-26-15-53','3491-10-08-12-20','3213-09-30-22-27','2946-04-06-17-40','3073-11-13-03-26','3007-09-29-22-15','2701-07-11-18-41','3321-06-15-15-40','2864-01-04-04-24','2591-11-25-02-37','2591-11-25-12-05','2509-08-28-22-49','2926-03-22-14-16','2539-01-17-04-44','2893-05-13-16-26','3025-07-30-21-26','3293-08-13-21-32','2605-01-29-22-44','3337-03-19-14-47','3157-09-07-14-06','3080-07-14-23-40','2862-05-10-13-41','2964-04-05-11-46','2757-07-15-21-26','3276-01-14-12-33'};
% exclude due to download issue, or multiple stays
%good_els = setxor(1:length(up.pt_id), [3,7,8,9,21,22,23,37, 39]);
%up.pt_id = up.pt_id(good_els);
%up.pt_data_segment = up.pt_data_segment(good_els);
up.req_data_time = 60*60;  % data taken from 1 hour into stay (secs)
up.req_data_duration = 8*60;  % 8-min segments of data extracted (secs)
up.dataset_name = 'bidmc';
up.csv.new_line = '\n';
up.units.signals = {'II', 'AVR', 'V',  'III', 'I',  'MCL1', 'MCL', 'PLETH', 'RESP', 'ABP',  'CVP',  'ART',   'P1',  'UAP',  'PAP'};
up.units.units   = {'mV', 'mV',  'mV', 'mV',  'mV', 'mV',   'mV',  'NU',    'pm'  , 'mmHg', 'mmHg', 'mmHg', 'mmHg', 'mmHg', 'mmHg'};
up.units.numerics    = {'HR',  'PULSE',  'SpO2', 'RESP', 'ABPSys', 'ABPDias', 'ABPMean', 'NBPSys', 'MBPDia', 'MBPMean'};
up.units.num_units   = {'bpm', 'bpm',    '%',    'pm',   'mmHg',   'mmHg',    'mmHg',    'mmHg',   'mmHg',   'mmHg'   };

% paths
up.paths.file_location = [up.paths.data_root, 'physionet.org', filesep, 'physiobank' ,filesep 'database', filesep];
up.paths.extracted_data_file = [up.paths.data_root, 'extracted_data.mat'];
up.paths.individual_files_folder = [up.paths.data_root, 'individual_files', filesep];

% url paths
up.paths.url.records = [up.paths.database_dir, '/RECORDS'];

% Add all functions within the directory
addpath(genpath(fileparts(mfilename('fullpath'))));

% check that all the required folders exist. If not, create them:
req_folders = {up.paths.data_root, up.paths.individual_files_folder};
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

function provide_details_of_script

licence_details = ['\n\n BIDMC_data_importer', ...
    '\n Copyright (C) 2017  King''s College London',...
    '\n This program comes with ABSOLUTELY NO WARRANTY', ... 
    '\n and is available under the GNU public license.', ...
    '\n For details see the accompanying LICENSE.txt file.\n\n'];

fprintf(licence_details)

end

function download_data(up)

fprintf('\n -- Downloading Data')

%% set current directory to that of wget:
cd(up.paths.data_root)

%% Download each record

for rec_no = 1 : length(up.pt_id)
    
    % specify record
    recording = [up.pt_id{rec_no}, '_', strrep(up.pt_data_segment{rec_no}, '-', '_')];
    
    % set save location
    save_folder = [up.paths.file_location, recording, filesep];
    if ~exist(save_folder, 'dir')
        mkdir(save_folder)
    end
    
    % download overall header file
    curr_file_name = [up.pt_id{rec_no}, '-', up.pt_data_segment{rec_no}, '.hea'];
    down_loc = [up.paths.database_dir, '/', up.pt_id{rec_no}, '/', curr_file_name];
    download_source.overall_header_file = down_loc;
    save_loc = [save_folder, curr_file_name];
    if ~exist(save_loc, 'file')
        outfilename = websave(save_loc,down_loc);
    end
    
    % find out which file contains the required data
    req_start_time = up.req_data_time;
    req_end_time = up.req_data_time + up.req_data_duration;
    
    % find out what individual files are in this record:
    fileID = fopen(save_loc);
    individual_files.name = {};
    individual_files.duration = []; individual_files.cum_duration_before = []; individual_files.cum_duration_after = [];
    line_no = 0;
    while ~feof(fileID)
        curr_text_line = fgetl(fileID); line_no = line_no+1;
        if line_no == 1
            temp = strfind(curr_text_line, ':');
            individual_files.start_time = curr_text_line(temp(1)-2:temp(1)+9);
        elseif ~isempty(strfind(curr_text_line, '_')) & isempty(strfind(curr_text_line, 'layout'))
            individual_files.name{end+1} = curr_text_line(1:strfind(curr_text_line, '_')+4);
            individual_files.duration(end+1) = str2double(curr_text_line(strfind(curr_text_line, '_')+6:end));
            individual_files.cum_duration_before(end+1) = sum(individual_files.duration(1:end-1));
            individual_files.cum_duration_after(end+1) = sum(individual_files.duration);
        elseif ~isempty(strfind(curr_text_line, '~')) & isempty(strfind(curr_text_line, 'layout'))
            individual_files.name{end+1} = curr_text_line(1:strfind(curr_text_line, '~'));
            individual_files.duration(end+1) = str2double(curr_text_line(strfind(curr_text_line, '~')+2:end));
            individual_files.cum_duration_before(end+1) = sum(individual_files.duration(1:end-1));
            individual_files.cum_duration_after(end+1) = sum(individual_files.duration);
        end
        clear curr_text_line
    end
    clear line_no
    fclose(fileID); clear fileID
    rel_file_el = find(individual_files.cum_duration_before <= (req_start_time*up.fs) & ...
        individual_files.cum_duration_after >= (req_end_time*up.fs)); clear req_start_time req_end_time
    rel_file_name = individual_files.name{rel_file_el};
    
    % download data header file
    curr_file_name = [rel_file_name, '.hea'];
    down_loc = [up.paths.database_dir, '/', up.pt_id{rec_no}, '/', curr_file_name];
    download_source.waves_header_file = down_loc;
    save_loc = [save_folder, curr_file_name];
    if ~exist(save_loc, 'file')
        outfilename = websave(save_loc,down_loc);
    end
    clear outfilename save_loc down_loc curr_file_name
    
    % also download data file
    curr_file_name = [rel_file_name, '.dat'];
    down_loc = [up.paths.database_dir, '/', up.pt_id{rec_no}, '/', curr_file_name];
    download_source.waves_data_file = down_loc;
    save_loc = [save_folder, curr_file_name];
    if ~exist(save_loc, 'file')
        outfilename = websave(save_loc,down_loc);
    end
    clear outfilename save_loc down_loc curr_file_name
    
    % download numerics header file
    curr_file_name = [up.pt_id{rec_no}, '-', up.pt_data_segment{rec_no}, 'n.hea'];
    down_loc = [up.paths.database_dir, '/', up.pt_id{rec_no}, '/', curr_file_name];
    save_loc = [save_folder, curr_file_name];
    if ~exist(save_loc, 'file')
        try
            outfilename = websave(save_loc,down_loc);
        catch
            current_ending = str2double(up.pt_data_segment{rec_no}(end-1:end));
            new_ending = current_ending-1;
            curr_file_name = [up.pt_id{rec_no}, '-', up.pt_data_segment{rec_no}(1:end-2), num2str(new_ending,'%02.f'), 'n.hea'];
            down_loc = [up.paths.database_dir, '/', up.pt_id{rec_no}, '/', curr_file_name];
            save_loc = [save_folder, curr_file_name];
            outfilename = websave(save_loc,down_loc);
        end
        download_source.numerics_header_file = down_loc;
    end
    clear outfilename save_loc down_loc curr_file_name
    
    % download numerics data file
    temp = strfind(rel_file_name, '_');
    curr_file_name = [rel_file_name(1:temp-1), 'n.dat']; clear temp
    down_loc = [up.paths.database_dir, '/', up.pt_id{rec_no}, '/', curr_file_name];
    download_source.numerics_data_file = down_loc;
    save_loc = [save_folder, curr_file_name];
    if ~exist(save_loc, 'file')
        outfilename = websave(save_loc,down_loc);
    end
    clear outfilename save_loc down_loc curr_file_name
    
    % save time duration at start of downloaded data file
    curr_file_name = [rel_file_name, '_Duration_Before.mat'];
    save_loc = [save_folder, curr_file_name];
    if ~exist(save_loc, 'file')
        duration_before = individual_files.cum_duration_before(rel_file_el);
        save(save_loc, 'duration_before', 'download_source');
    end
    clear waveforms_start_time duration_before save_loc curr_file_name rel_file_name rel_file_el
    
end

end

function extract_from_mimic_ii(up)
%
% Based on RDSAMP, part of the WFDB Toolbox. Further details are available
% at:
%       http://physionet.org/physiotools/matlab/wfdb-app-matlab/
%
% Reads MIMIC II data and outputs it as a structure of .t and .v values.

fprintf('\n -- Extracting Data from MIMIC II files')

for pt_stay_no = 1 : length(up.pt_id)
    
    % set current dir
    dir_path = [up.paths.file_location, up.pt_id{pt_stay_no}, '_', strrep(up.pt_data_segment{pt_stay_no}, '-', '_'), filesep];
    cd(dir_path)
    
    % setup java
    % persistent javaWfdbExec config
    [javaWfdbExec, config] = deal([]);
    if(isempty(javaWfdbExec))
        [javaWfdbExec,config]=getWfdbClass('rdsamp');
    end
    
    % Remove file extension
    temp = dir([dir_path, '*.dat']);
    names = extractfield(temp, 'name');
    names = names(~cellfun(@isempty, strfind(names, '_'))); names = names{1,1};
    record_str = names(1:end-4); clear temp names
    
    % Extract file information (including which signals are present)
    siginfo = extract_file_info(record_str);
    
    % Load duration at start of file
    load([dir_path, record_str, '_Duration_Before.mat'])
    
    % create wfdb argument
    N0 = (siginfo.fs * up.req_data_time) - duration_before;
    wfdb_argument={'-r',record_str,'-Ps','-f',['s' num2str(N0)]}; clear record_str
    
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
    
    % store the signals
    for wave_no = 1 : length(siginfo.sig_names)
        eval(['waves.' siginfo.sig_names{wave_no} '.t = extracted_data(:, 1)-(N0/125);']);
        eval(['waves.' siginfo.sig_names{wave_no} '.v = extracted_data(:, wave_no+1);']);
    end
    clear wave_no N N0 current_waves extracted_data Fs ListCapacity signalList javaWfdbExec config
    
    
    %% Extract corresponding numerics
        
    % Extract file information (including which signals are present)
    header_filepath = dir([dir_path, '*n.hea']);
    if ~isempty(header_filepath)
        header_filepath = header_filepath.name(1:end-4);
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
        start_time_diff = (datenum(['00/01/0000 ' numinfo.start_time], 'dd/mm/yyyy HH:MM:SS') - datenum(['00/01/0000 ' siginfo.start_time], 'dd/mm/yyyy HH:MM:SS')) *24*60*60;  % in secs
        N0 = (numerics.fs * up.req_data_time) - (duration_before/125) - start_time_diff;
        clear duration_before
        N = N0 +(numerics.fs * up.req_data_duration);
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
    
    overall_header_file = [dir_path, up.pt_id{pt_stay_no}, '-', up.pt_data_segment{pt_stay_no}, '.hea'];
    fileID = fopen(overall_header_file);
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
    
    %% Add annotations
    ann = add_annotations(pt_stay_no);
    
    % Load from file
    ref_files_folder = '/Users/petercharlton/Downloads/data_from_marco/';
    ref_file_name = [up.pt_id{pt_stay_no}, '_', strrep(up.pt_data_segment{pt_stay_no}, '-', '_'), '.mat'];
    load([ref_files_folder, ref_file_name], 'ann'); clear ref_file_name ref_files_folder
    
    %% Store data for this patient stay in the output
    data{pt_stay_no}.loc = fix.loc;
    data{pt_stay_no}.age = fix.age;
    data{pt_stay_no}.sex = fix.sex;
    data{pt_stay_no}.ID = up.pt_id{pt_stay_no};
    data{pt_stay_no}.data_segment = up.pt_data_segment{pt_stay_no};
    data{pt_stay_no}.waves = waves;
    data{pt_stay_no}.ann = ann;
    data{pt_stay_no}.download_source = download_source;
    clear waves numerics fix ann download_source
    
end

%% Save extracted data
save(up.paths.extracted_data_file, 'data')

end

function save_data_to_individual_files(up)

fprintf('\n -- Saving data to individual files')

%% Load extracted data
load(up.paths.extracted_data_file, 'data')
all_data = data; clear data

%% Save to individual files
for file_no = 1 : length(all_data)
    
    % specify filepath to save this patient stay data in
    filepath = [up.paths.individual_files_folder, 'Ind_data_' all_data{file_no}.ID, '_', strrep(all_data{file_no}.data_segment, '-', '_')];
    
    % create fixed variables to be saved
    age = all_data{file_no}.age;
    gender = all_data{file_no}.sex;
    location = all_data{file_no}.loc;
    ventMode = '';
    download_source = all_data{file_no}.download_source;
    
    % create variable relating to signals
    signals = fieldnames(all_data{file_no}.waves);
    signals = signals(~strcmp(signals, 'fs'));
    fs = repmat(all_data{file_no}.waves.fs, [1, length(signals)]);
    header = signals(:)';
    data = [];
    for wave_type_no = 1 : length(signals)
        eval(['data(:,end+1) = all_data{file_no}.waves.' signals{wave_type_no} '.v;'])
    end
    clear wave_type_no signals
    
    % copy across annotations
    ann = all_data{file_no}.ann;
    
    % save variables to file
    save(filepath, 'age', 'data', 'fs', 'gender', 'header', 'location', 'ventMode', 'ann', 'download_source');
    
    % save numerics if they exist
    if sum(strcmp(fieldnames(all_data{file_no}), 'numerics'))
        numerics = all_data{file_no}.numerics;
        save(filepath, 'numerics', '-append')
    end
    clear age data numerics fs gender header location ventMode filepath ann
    
end
clear file_no


end

function check_individual_files(up)

% setup
ref_files_folder = '/Users/petercharlton/Downloads/data_from_marco/';

fprintf('\n -- Checking Extracted Data')

% create list of files
temp = dir([up.paths.individual_files_folder, '*.mat']);
created_files = extractfield(temp, 'name'); clear temp
[stats.missing, stats.different] = deal(nan(length(created_files),1));

for file_no = 1 : length(created_files)
    curr_file = created_files{file_no};
    created_path = [up.paths.individual_files_folder, curr_file];
    ref_path = [ref_files_folder, curr_file(10:end)];
    clear curr_file
    
    % load the two data files
    created_data = load(created_path);
    ref_data = load(ref_path);
    clear ref_path created_path
    
    % setup storage
    [temp.missing_vars, temp.different_vars] = deal({});
    
    created_variables = fieldnames(created_data);
    
    % check each fixed variable in turn:
    ref_variables = fieldnames(ref_data);
    ref_variables = ref_variables(~strcmp(ref_variables, 'data'));
    ref_variables = ref_variables(~strcmp(ref_variables, 'fs'));
    ref_variables = ref_variables(~strcmp(ref_variables, 'header')); % up.extraction.rel_sigs
    for var_no = 1 : length(ref_variables)
        curr_var = ref_variables{var_no};
        if ~sum(strcmp(created_variables, curr_var))
            temp.missing_vars(end+1) = curr_var;
            continue
        end
        eval(['ref_var = ref_data.' curr_var ';']);
        eval(['created_var = created_data.' curr_var ';']);
        if strcmp(curr_var, 'data') || strcmp(curr_var, 'header')
            [~, rel_ref_order] = sort(ref_data.header);
            [~, rel_created_order] = sort(created_data.header);
            ref_var = ref_var(:,rel_ref_order);
            created_var = created_var(:,rel_created_order);
            clear rel_created_order rel_ref_order
        end
        if isequal(created_var, ref_var) || (isnan(ref_var) && (isnumeric(created_var) && isnan(created_var)))
            continue
        else
            temp.different_vars(end+1) = {curr_var};
        end
    end
    clear var_no curr_var ref_var created_var
    
    % check each signal-related variable in turn
    for curr_var = {'fs', 'header', 'data'}
        switch curr_var{1,1}
            case 'fs'
                created_var = created_data.fs(1);
                ref_var = ref_data.fs(1);
            case 'header'
                [created_var,c1,c2] = intersect(created_data.header, up.extraction.rel_sigs);
                [ref_var,r1,r2] = intersect(ref_data.header, up.extraction.rel_sigs);
            case 'data'
                created_var = created_data.data(:,c1);
                ref_var = ref_data.data(:,r1);
        end
        if isequal(created_var, ref_var)
            continue
        else
            temp.different_vars(end+1) = curr_var;
        end
    end
    
    % calculate stats
    stats.missing(file_no) = length(temp.missing_vars);
    stats.different(file_no) = length(temp.different_vars);
    if stats.different(file_no) > 0
        subplot(3,1,1)
        plot(created_var(1:1000,1), '.-'), hold on, plot(ref_var(1:1000,1), 'o-')
        subplot(3,1,2)
        plot(created_var(1:1000,2), '.-'), hold on, plot(ref_var(1:1000,2), 'o-')
        subplot(3,1,3)
        plot(created_var(1:1000,3), '.-'), hold on, plot(ref_var(1:1000,3), 'o-')
        title(temp.different_vars)
        temp.different_vars
        waitfor(gcf)
    end
    clear temp
end
clear file_no

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
headers.start_time = temp.top_line{rel_line}; clear rel_line

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

function ann = add_annotations(subj_no)

switch subj_no
    case 1
        ann.resp.breaths.x1 = [189,476,737,1037,1324,1646,2054,2350,2628,2924,3220,3528,3848,4191,4517,4883,5226,5565,5922,6283,6626,7004,7409,7776,8128,8507,8841,9193,9533,9907,10250,10567,10876,11202,11526,11830,12143,12483,12835,13170,13517,13861,14217,14583,14930,15272,15641,16024,16393,16776,17146,17511,17854,18189,18515,18857,19196,19526,19878,20235,20617,21000,21330,21691,22052,22417,22785,23159,23498,23841,24172,24498,24867,25220,25576,25915,26287,26648,26996,27343,27722,28087,28461,28826,29183,29548,29909,30237,30602,30959,31341,31685,32102,32446,32807,33176,33541,33909,34274,34639,35009,35387,35709,36087,36452,36783,37152,37513,37846,38224,38559,38954,39324,39672,40046,40420,40802,41167,41526,41870,42230,42587,42943,43287,43704,44139,44483,44843,45189,45567,45889,46267,46663,46998,47354,47707,48063,48428,48809,49139,49500,49865,50200,50543,50904,51243,51626,51991,52335,52672,53007,53320,53672,54024,54398,54763,55098,55467,55807,56167,56509,56870,57252,57626,57974,58352,58704,59078,59426,59809];
        ann.resp.breaths.x2 = [187,472,747,1042,1322,1649,2061,2346,2625,2910,3222,3522,3842,4180,4513,4892,5230,5557,5937,6296,6634,7035,7436,7788,8131,8526,8864,9202,9540,9925,10263,10590,10906,11223,11548,11844,12150,12487,12841,13189,13527,13865,14229,14582,14941,15272,15652,16021,16391,16797,17140,17536,17873,18206,18528,18858,19191,19528,19892,20246,20615,21001,21333,21713,22051,22420,22793,23157,23505,23854,24170,24518,24861,25231,25590,25922,26268,26659,26997,27340,27730,28094,28464,28828,29166,29551,29894,30245,30620,30958,31338,31702,32092,32441,32831,33185,33570,33911,34286,34655,35030,35378,35726,36090,36454,36792,37146,37513,37856,38210,38584,38954,39344,39677,40057,40426,40806,41191,41532,41870,42224,42598,42936,43290,43696,44160,44498,44841,45203,45562,45910,46264,46649,46992,47367,47710,48058,48438,48784,49133,49486,49877,50209,50557,50900,51280,51629,51972,52357,52677,53004,53331,53674,54017,54392,54756,55109,55479,55822,56165,56511,56849,57245,57651,57994,58332,58728,59071,59414,59809];
    case 2
        ann.resp.breaths.x1 = [311,767,1246,1741,2228,2689,3146,3641,4122,4583,5039,5513,5987,6443,6917,7378,7850,8328,8811,9220,9763,10267,10798,11317,11765,12287,12817,13291,13813,14317,14813,15293,15789,16276,16802,17280,17741,18207,18654,19135,19622,20100,20778,21309,21843,22348,22815,23315,23824,24346,24798,25250,25750,26278,26774,27270,27726,28213,28700,29170,29674,30198,30702,31167,31633,32150,32633,33089,33667,34152,34700,35209,35717,36239,36717,37222,37689,38185,38667,39154,39624,40076,40341,40576,40959,41400,41883,42365,42887,43426,43891,44426,44900,45324,45998,46372,46841,47311,47815,48402,48952,49509,50039,50548,51083,51657,52161,52754,53224,53720,54193,54650,55128,55589,56098,56496,57035,57548,58087,58613,59091,59552];
        ann.resp.breaths.x2 = [319,784,1264,1755,2251,2689,3174,3665,4127,4586,5051,5526,6001,6465,6919,7394,7888,8342,8817,9234,9761,10284,10801,11311,11780,12297,12825,13311,13796,14334,14830,15303,15805,16290,16818,17282,17773,18201,18681,19138,19634,20135,20811,21333,21850,22362,22830,23326,23843,24360,24814,25289,25764,26268,26780,27271,27725,28226,28712,29197,29688,30187,30710,31195,31644,32150,32641,33079,33665,34154,34697,35199,35726,36254,36734,37230,37703,38204,38679,39165,39624,40099,40326,40590,40959,41422,41907,42398,42899,43416,43896,44408,44894,45293,46000,46375,46834,47319,47789,48354,48974,49497,50009,50552,51096,51660,52177,52756,53241,53727,54181,54645,55099,55595,56102,56501,57023,57556,58089,58617,59108,59556];
    case 3
        ann.resp.breaths.x1 = [293,746,1198,1646,2085,2511,2880,3315,3783,4222,4670,5074,5474,5891,6339,6748,7152,7611,8041,8415,8850,9285,9689,10124,10567,10976,11400,11857,12274,12691,13196,13491,13874,14296,14726,15163,15554,16007,16420,16824,17207,17576,17985,18385,19148,19522,19922,20326,20739,21157,21539,21957,22413,22824,23237,23685,24098,24559,24980,25437,25828,26246,26652,27113,27509,27891,28339,28700,29122,29604,30002,30433,30820,31237,31685,32076,32580,32998,33420,33813,34265,34635,35074,35457,35865,36309,36730,37161,37646,38089,38524,38920,39233,39693,40115,40572,40980,41378,41843,42230,42635,43057,43457,43891,44287,44817,45241,45663,46050,46472,46902,47324,47724,48098,48502,48900,49335,49687,50087,50522,50922,51426,51891,52330,52793,53233,53698,54163,54585,55028,55476,55876,56278,56691,57143,57609,58048,58509,58952,59422,59887];
        ann.resp.breaths.x2 = [314,747,1195,1612,2092,2509,2889,3332,3768,4243,4676,5067,5462,5916,6354,6750,7151,7640,8030,8426,8827,9276,9703,10099,10563,10970,11400,11838,12234,12704,13173,13490,13881,14345,14751,15145,15541,16005,16412,16802,17193,17588,17968,18364,18753,19148,19518,19924,20320,20732,21133,21534,21940,22410,22819,23215,23690,24080,24571,24972,25368,25848,26239,26638,27102,27509,27910,28305,28712,29139,29593,30008,30409,30789,31206,31686,32092,32578,32984,33380,33784,34275,34655,35067,35462,35858,36265,36739,37151,37635,38120,38505,38912,39265,39693,40099,40574,40980,41369,41838,42239,42625,43031,43432,43918,44292,44799,45251,45662,46063,46475,46913,47346,47726,48121,48522,48911,49312,49692,50082,50526,50937,51428,51882,52336,52814,53215,53701,54181,54582,55041,55474,55875,56268,56664,57150,57635,58047,58532,58939,59419,59920];
    case 4
        ann.resp.breaths.x1 = [202,663,1020,1341,1720,2211,2663,3141,3567,4057,4522,4909,5378,5878,6296,6765,7235,7659,8098,8602,8972,9463,9967,10450,10876,11296,11774,12161,12670,13117,13600,13991,14378,14870,15346,15759,16276,16741,17128,17528,17972,18385,18804,19252,19722,20117,20565,21061,21530,21935,22313,22802,23276,23676,24150,24550,25093,25624,26085,26383,26874,27309,27752,28235,28726,29139,29561,30011,30485,30959,31380,31989,32493,32998,33424,33887,34304,34730,35191,35674,36135,36548,37048,37528,38007,38489,38880,39359,39850,40163,40507,40920,41246,41635,42043,42522,43030,43426,43874,44339,44748,45180,45672,46154,46763,47259,47711,48163,48680,49170,49565,50022,50457,50909,51343,51783,52287,52746,53246,53711,54241,54676,55089,55611,56054,56509,56957,57409,57957,58474,58904,59370,59861];
        ann.resp.breaths.x2 = [198,662,1021,1338,1718,2214,2673,3153,3570,4080,4528,4929,5394,5885,6296,6771,7235,7666,8094,8606,8964,9466,9962,10447,10885,11284,11770,12176,12688,13131,13606,13986,14382,14862,15340,15757,16285,16744,17140,17520,17968,18406,18795,19254,19724,20125,20568,21064,21528,21924,22315,22814,23284,23669,24165,24550,25083,25627,26075,26395,26865,27303,27751,28232,28728,29139,29556,30008,30488,30958,31380,32003,32504,32995,33433,33922,34275,34734,35188,35679,36133,36560,37046,37540,38004,38484,38896,39344,39856,40162,40500,40917,41239,41643,42065,42546,43020,43437,43891,44345,44751,45203,45683,46169,46760,47261,47704,48164,48686,49185,49560,50019,50468,50911,51349,51787,52288,52745,53252,53722,54244,54682,55104,55632,56049,56506,56949,57408,57957,58474,58923,59377,59867];
    case 5
        ann.resp.breaths.x1 = [67,1285,2572,3839,5035,6109,7615,8902,10167,11335,12291,13574,14904,16059,17520,18693,20483,21791,22841,23772,24815,26033,27435,28926,30093,31424,32480,33433,34561,35717,37161,38837,40450,41613,42922,44174,45293,46980,48185,49543,50578,51722,53076,54524,55720,56857,57943,58991];
        ann.resp.breaths.x2 = [55,1264,2578,3842,5051,6117,7645,8880,10173,11284,12276,13564,14925,16053,17483,18644,20468,21771,22851,23759,24809,26044,27429,28923,30055,31417,32467,33427,34555,35726,37167,38796,40453,41585,42904,44171,45277,46960,48158,49539,50594,51703,53057,54518,55706,56865,57915,58970];
    case 6
        ann.resp.breaths.x1 = [380,741,1133,1485,1893,2250,2624,3007,3380,3761,4126,4496,4865,5243,5630,6009,6383,6770,7130,7498,7889,8228,8611,8993,9367,9741,10111,10515,10841,11163,11613,11991,12339,12752,13135,13500,13865,14257,14622,15000,15376,15754,16115,16493,16867,17263,17611,17998,18359,18746,19104,19487,19874,20239,20630,20983,21378,21752,22117,22504,22863,23250,23620,23998,24363,24750,25128,25493,25876,26246,26613,26991,27383,27730,28126,28483,28883,29261,29604,30011,30376,30737,31115,31515,31846,32246,32624,32998,33372,33759,34117,34491,34878,35252,35626,36000,36374,36739,37104,37504,37854,38237,38628,38989,39363,39746,40115,40493,40867,41259,41617,41996,42357,42761,43130,43509,43887,44248,44609,45000,45376,45733,46120,46480,46841,47254,47615,47976,48367,48746,49087,49087,49509,49883,50230,50596,50996,51374,51743,52122,52487,52846,53237,53620,53998,54376,54746,55115,55498,55885,56241,56613,56983,57370,57739,58117,58491,58874,59257,59609,59996];
        ann.resp.breaths.x2 = [372,747,1121,1491,1865,2240,2615,2995,3375,3749,4127,4491,4866,5241,5615,5990,6370,6750,7130,7494,7872,8236,8621,8991,9371,9735,10109,10500,10859,11228,11617,11997,12303,12741,13115,13490,13865,14245,14625,14984,15372,15741,16127,16485,16865,17235,17620,17989,18369,18739,19122,19491,19861,20235,20615,20985,21365,21745,22114,22508,22861,23247,23611,23991,24365,24740,25115,25484,25859,26239,26606,26986,27366,27735,28115,28485,28865,29245,29614,29999,30367,30741,31121,31491,31844,32235,32609,32984,33369,33739,34127,34497,34871,35246,35600,36001,36365,36750,37114,37494,37861,38236,38611,38991,39371,39735,40109,40495,40864,41244,41617,41991,42361,42746,43121,43495,43865,44224,44614,44978,45356,45747,46116,46491,46860,47245,47615,47979,48369,48749,49117,49491,49866,50235,50605,50985,51370,51745,52119,52489,52872,53231,53616,53980,54360,54740,55120,55489,55864,56244,56617,56991,57361,57730,58110,58485,58859,59229,59614,59989];
    case 7
        ann.resp.breaths.x1 = [180,563,928,1315,1672,2050,2446,2785,3163,3567,3926,4309,4678,5052,5417,5813,6170,6583,6909,7309,7676,8054,8428,8807,9189,9541,9941,10298,10676,11063,11430,11796,12174,12539,12957,13317,13691,14048,14448,14800,15180,15537,15933,16298,16680,17050,17437,17789,18185,18528,18935,19304,19661,20048,20430,20809,21170,21565,21909,22309,22680,23033,23476,23820,24163,24537,24924,25320,25685,26041,26417,26791,27183,27552,27930,28326,28678,29057,29470,29800,30189,30554,30941,31293,31672,32033,32420,32807,33167,33533,33939,34304,34687,35061,35413,35796,36183,36543,36935,37304,37667,38076,38446,38785,39176,39537,39924,40302,40689,41020,41435,41796,42178,42530,42939,43313,43739,44074,44413,44796,45167,45563,45924,46307,46667,47037,47446,47889,48167,48541,48909,49296,49678,50035,50430,50804,51178,51548,51917,52309,52672,53059,53424,53789,54167,54541,54920,55293,55685,56059,56439,56809,57170,57552,57939,58313,58661,59048,59417,59804];
        ann.resp.breaths.x2 = [182,573,937,1306,1675,2050,2425,2805,3179,3559,3922,4307,4676,5061,5431,5800,6175,6555,6935,7309,7677,8051,8426,8801,9175,9550,9925,10299,10690,11049,11432,11796,12181,12546,12941,13305,13675,14044,14435,14804,15177,15551,15921,16301,16681,17045,17420,17799,18174,18554,18927,19301,19687,20046,20436,20800,21180,21549,21929,22299,22672,23046,23442,23801,24181,24555,24930,25299,25679,26049,26432,26801,27176,27546,27931,28305,28675,29049,29429,29809,30177,30551,30926,31296,31675,32045,32435,32815,33179,33544,33927,34301,34676,35056,35425,35805,36180,36549,36935,37299,37666,38046,38442,38801,39175,39561,39920,40305,40679,41033,41427,41801,42171,42551,42925,43300,43717,44055,44429,44804,45172,45546,45916,46322,46670,47045,47425,47868,48174,48554,48916,49312,49676,50051,50425,50800,51175,51549,51919,52299,52661,53051,53416,53801,54175,54545,54930,55294,55664,56049,56432,56801,57181,57540,57925,58300,58675,59039,59419,59799];
    case 8
        ann.resp.breaths.x1 = [254,511,798,1215,1567,1885,2202,2511,2828,3193,3524,3917,4313,4678,5039,5491,5813,6170,6578,6965,7309,7633,7985,8350,8720,9111,9524,9859,10176,10593,10933,11252,11583,11961,12326,12648,12978,13352,13674,14013,14322,14661,14943,15328,15707,16067,16459,16833,17215,17593,17876,18259,18533,18922,19330,19704,20074,20430,20730,21091,21561,21900,22261,22611,22976,23298,23628,23959,24341,24737,25111,25480,25859,26172,26430,26743,27143,27517,27874,28170,28539,28917,29239,29613,29935,30320,30698,31046,31389,31720,32076,32446,32741,33098,33485,33822,34187,34509,34835,35209,35509,35843,36222,36574,36943,37317,37663,38076,38454,38833,39207,39628,39959,40272,40585,40911,41265,41648,42004,42400,42696,43017,43365,43678,44048,44365,44735,45080,45437,45798,46163,46489,46837,47176,47533,47841,48267,48576,48952,49335,49761,50096,50522,50913,51270,51600,51904,52283,52663,53072,53402,53754,54046,54376,54763,55154,55515,55863,56172,56491,56791,57152,57483,57804,58239];
        ann.resp.breaths.x2 = [256,509,794,1206,1575,1887,2235,2520,2826,3206,3517,3895,4307,4671,5035,5494,5816,6191,6581,6966,7304,7650,7993,8342,8727,9096,9503,9846,10183,10595,10933,11258,11590,11949,12313,12656,12999,13348,13675,14013,14313,14635,14957,15340,15710,16069,16470,16813,17219,17588,17884,18269,18570,18927,19323,19697,20061,20415,20753,21117,21549,21898,22251,22613,22951,23299,23621,23954,24355,24719,25104,25474,25854,26165,26427,26770,27134,27524,27883,28184,28548,28928,29250,29614,29957,30319,30694,31037,31380,31734,32071,32441,32736,33090,33480,33827,34175,34507,34834,35209,35510,35842,36217,36586,36956,37330,37656,38057,38463,38838,39202,39608,39946,40278,40574,40912,41258,41638,42028,42398,42704,43010,43358,43691,44039,44366,44735,45066,45441,45805,46179,46485,46823,47161,47509,47842,48269,48596,48953,49338,49745,50098,50547,50906,51265,51613,51898,52294,52672,53078,53405,53748,54054,54376,54756,55131,55521,55864,56170,56495,56791,57144,57482,57809,58258,58527,58907,59261,59651,59989];
    case 9 
        ann.resp.breaths.x1 = [320,685,1080,1450,1802,2193,2559,2946,3333,3685,4061,4426,4822,5200,5552,5917,6309,6700,7052,7417,7820,8198,8554,8933,9320,9702,10067,10437,10815,11193,11574,11943,12326,12683,13078,13443,13813,14187,14565,14948,15315,15680,16072,16446,16815,17189,17567,17946,18307,18698,19070,19448,19817,20183,20578,20930,21296,21691,22074,22426,22811,23228,23554,23941,24315,24685,25072,25441,25828,26189,26565,26939,27322,27687,28061,28435,28809,29191,29565,29939,30307,30685,31054,31411,31811,32211,32524,32915,33289,33676,34052,34448,34813,35178,35552,35939,36317,36687,37052,37461,37811,38193,38572,38941,39328,39702,40063,40433,40802,41185,41557,41961,42304,42678,43048,43435,43813,44191,44570,44935,45315,45689,46054,46441,46802,47185,47572,47924,48315,48672,49078,49435,49796,50191,50557,50939,51304,51687,52057,52426,52815,53211,53537,53928,54302,54685,55050,55446,55802,56176,56561,56948,57317,57696,58061,58430,58804,59170,59574,59952];
        ann.resp.breaths.x2 = [325,694,1074,1443,1823,2198,2567,2947,3317,3702,4075,4449,4819,5188,5568,5948,6323,6703,7072,7447,7819,8199,8563,8949,9323,9703,10073,10447,10822,11191,11569,11954,12329,12693,13073,13443,13817,14197,14577,14947,15319,15694,16063,16443,16828,17193,17573,17947,18317,18691,19075,19454,19819,20199,20563,20937,21323,21697,22072,22447,22814,23189,23574,23949,24323,24693,25062,25442,25822,26191,26564,26939,27319,27693,28063,28443,28812,29192,29572,29941,30314,30683,31069,31433,31823,32203,32567,32947,33317,33686,34069,34439,34808,35188,35568,35943,36323,36681,37067,37447,37819,38194,38563,38933,39323,39698,40062,40442,40811,41197,41569,41939,42319,42688,43057,43427,43812,44197,44567,44947,44340,45309,45689,46069,46443,46818,47187,47562,47942,48327,48697,49064,49444,49808,50188,50573,50937,51317,51692,52061,52441,52819,53189,53569,53943,54307,54693,55073,55442,55817,56191,56559,56939,57313,57693,58068,58437,58817,59197,59567,59936];
    case 10
        ann.resp.breaths.x1 = [293,663,1015,1346,1672,1998,2328,2667,3020,3420,3857,4322,4713,5148,5535,5913,6326,6726,7117,7500,7872,8263,8637,9059,9493,9920,10354,10715,11115,11509,11900,12322,12739,13065,13443,13852,14313,14726,15163,15567,15924,16337,16720,17102,17485,17898,18302,18702,19109,19548,19917,20313,20878,21265,21648,22048,22430,22820,23202,23563,23985,24380,24767,25128,25507,25907,26296,26730,27152,27535,27957,28352,28752,29135,29574,30000,30389,30772,31485,31950,32363,32820,33228,33611,34057,34478,34887,35274,35730,36174,36535,36939,37248,37841,38250,38607,38933,39372,39659,40020,40372,40789,41024,41448,41809,42183,42578,42987,43378,43761,44243,44548,45080,45515,45933,46380,46772,47241,47620,48107,48515,48961,49352,49761,50178,50574,50952,51374,51804,52243,52646,53059,53476,54163,54598,54998,55454,55841,56241,56683,57283,57674,58122,58487,58822,59596];
        ann.resp.breaths.x2 = [293,689,1042,1322,1644,2018,2330,2668,3032,3427,3858,4312,4724,5146,5536,5922,6328,6718,7114,7513,7872,8263,8653,9054,9492,9925,10326,10722,11107,11511,11907,12319,12720,13078,13464,13859,14303,14725,15156,15562,15921,16343,16712,17098,17493,17889,18296,18712,19096,19549,19908,20294,20563,20879,21270,21655,22030,22425,22819,23189,23579,23980,24376,24766,25168,25505,25896,26300,26733,27129,27546,27947,28358,28749,29134,29588,30013,30383,30778,31153,31485,31934,32372,32815,33222,33623,34043,34481,34866,35230,35721,36154,36534,36703,36898,37235,37840,38236,38584,38938,39371,39650,40009,40400,40785,41012,41416,41791,42171,42556,42947,43363,43765,44203,44598,45050,45493,45931,46369,46770,47172,47594,48100,48507,48932,49354,49724,50162,50557,50974,51375,51792,52214,52640,53046,53463,53848,54154,54582,54993,55437,55838,56239,56648,56949,57303,57699,58094,58480,58807,59150,59582];
    case 11
        ann.resp.breaths.x1 = [324,989,1615,1954,2498,3020,3546,4048,4513,5017,5474,5917,6487,7009,7820,8367,8850,9350,9907,10311,10828,11335,11796,12304,12804,13396,13922,14443,14948,15254,15585,16093,16372,17311,17833,18376,18996,19378,19987,20257,20783,21183,21670,22170,22733,23259,23850,24376,24941,25502,25933,26461,26917,27396,27743,28517,28904,29743,30254,30776,31285,31824,32411,32920,33472,33948,34422,34943,35461,35970,36535,37061,37585,38093,38641,39167,39650,40076,40767,41600,42048,42526,42983,43517,44013,44530,45028,45467,45959,46398,46702,46989,47454,47911,48428,48970,49530,50026,50561,51117,51600,52143,52667,53207,53750,54263,54885,55411,55915,56396,56861,57426,58039,58678,59113,59661];
        ann.resp.breaths.x2 = [309,1000,1675,1987,2530,2995,3311,3501,3980,4486,5135,5462,5932,6454,6940,7809,8373,8843,9328,9935,10273,10774,11279,11796,12271,12788,13406,13896,14424,14920,15266,15578,16069,16380,17293,17799,18348,19011,19349,19998,20199,20747,21164,21681,22141,22745,23273,23832,24365,24930,25521,25901,26400,26875,27398,27757,28506,28891,29693,30224,30799,31269,31839,32435,32900,33496,33911,34391,34956,35462,36016,36518,37082,37566,38099,38637,39165,39640,40083,40737,41128,41559,42028,42551,42962,43516,44013,44530,45045,45462,45974,46243,46406,46718,47008,47435,47900,48121,48427,48958,49513,50072,50578,51138,51586,52130,52629,53168,53716,54255,54883,55373,55922,56406,56881,57456,58063,58675,59203,59656];
    case 12
        ann.resp.breaths.x1 = [124,537,915,1354,1733,2176,2563,2941,3363,3783,4200,4587,5030,5465,5891,6348,6748,7157,7576,7989,8367,8754,9141,9511,9911,10320,10720,11107,11552,11930,12313,12704,13117,13504,13891,14300,14678,15072,15498,15941,16337,16728,17128,17572,17972,18393,18796,19209,19587,20013,20391,20800,21183,21639,22035,22448,22859,23298,23707,24115,24698,25098,25537,25976,26370,26926,27348,27752,28161,28591,29039,29452,29952,30498,30915,31424,31828,32228,32628,33002,33415,33809,34213,34570,34930,35283,35661,36030,36378,36730,37130,37487,37828,38185,38524,38902,39324,39672,40041,40415,40767,41124,41474,41852,42209,42552,42935,43296,43687,44065,44404,44787,45159,45528,45898,46237,46615,46972,47328,47680,48063,48411,48791,49135];
        ann.resp.breaths.x2 = [119,509,926,1348,1755,2156,2567,2958,3369,3768,4212,4592,5014,5457,5922,6333,6745,7151,7566,7978,8389,8743,9139,9513,9909,10326,10716,11128,11553,11912,12297,12693,13115,13501,13875,14287,14677,15071,15483,15926,16327,16734,17135,17573,17974,18391,18779,19185,19592,20003,20415,20795,21175,21618,22061,22468,22888,23294,23722,24117,24708,25094,25547,25970,26400,26918,27340,27751,28173,28596,28997,29435,29936,30493,30905,31401,31813,32224,32625,33005,33417,33805,34212,34576,34924,35294,35658,36022,36381,36745,37109,37478,37846,38178,38526,38917,39313,39672,40041,40400,40764,41117,41474,41844,42187,42556,42941,43290,43675,44055,44398,44772,45150,45551,45879,46248,46617,46976,47330,47689,48053,48417,48774,49127,49497,49856,50235,50610,50980,51344,51761,52172,52513,52914,53273,53648,54012,54408,54782,55220,55606,56007,56374,56754,57113,57503,57920,58327,58733,59134,59540,59973];
    case 13
        ann.resp.breaths.x1 = [120,537,898,1272,1659,2063,2459,2841,3228,3628,4000,4413,4809,5196,5630,6026,6383,6761,7204,7607,7998,8428,8850,9233,9602,10015,10420,10780,11124,11517,11935,12300,12743,13126,13470,13857,14191,14661,15072,15454,15780,16172,16533,16898,17289,17715,18054,18463,18848,19252,19630,20030,20443,20817,21213,21609,22009,22387,22750,23154,23493,23841,24272,24598,24993,25450,25863,26163,26630,26987,27335,27709,28074,28517,28857,29261,29643,30072,30446,30824,31272,31659,31980,32385,32763,33150,33524,33904,34248,34648,35022,35339,35752,36230,36622,36930,37300,37663,38007,38272,38663,38976,39315,39685,39993,40341,40698,41037,41335,41670,42017,42387,42722,43026,43352,43700,44048,44383,44704,45098,45476,45837,46189,46528,46859,47211,47550,47902,48311,48646,49030,49387,49770,50100,50474,50826,51139,51557,51870,52209,52924,53254,53620,53972,54385,54759,55128,55441,55820,56137,56557,56965,57326,57678,58052,58422,58765,59096,59478,59843];
        ann.resp.breaths.x2 = [61,277,515,858,1069,1280,1633,2077,2430,2858,3206,3596,3980,4396,4776,5204,5652,5985,6328,6776,7125,7468,7566,8015,8363,8796,9260,9587,9983,10368,10722,11107,11495,11896,12292,12730,13063,13406,13854,14208,14656,14994,15108,15425,15768,16127,16475,16929,17272,17710,18047,18427,18842,19191,19629,19977,20425,20763,21206,21544,22009,22341,22682,23125,23474,23811,24244,24577,25020,25379,25822,26170,26611,26949,27297,27683,28068,28516,28865,29187,29641,29968,30425,30757,31216,31544,31997,32335,32720,33111,33459,33906,34233,34592,35019,35352,35716,36243,36597,36929,37267,37598,37951,38278,38632,38959,39307,39645,39988,40315,40653,40991,41212,41675,41991,42334,42672,42999,43321,43670,44013,44350,44672,45010,45372,45815,46142,46485,46818,47156,47483,47852,48264,48591,48943,49386,49739,50056,50404,50742,51090,51528,51882,52204,52534,52872,53226,53553,53896,54334,54687,55120,55468,55806,56139,56480];
    case 14
        ann.resp.breaths.x1 = [267,841,1324,1937,2485,3154,3620,4252,4539,4865,5196,5465,6000,6535,7174,7728,8267,8854,9472,10098,10641,11326,11861,12491,12978,13574,14174,14709,15233,15898,16528,17024,17572,18167,18774,19317,19835,20417,21087,21609,22217,22854,23459,23946,24250,24811,25350,25980,26443,27035,27587,28165,28826,29400,29991,30467,31085,31637,32198,32750,33454,34013,34570,35148,35678,36257,36809,37326,37650,37928,38302,38563,39311,39846,40367,40989,41478,42100,42678,43274,43861,44370,44896,45520,46107,46654,47072,47324,47963,48637,49152,49739,50243,50783,51448,52009,52580,53050,53728,54367,54963,55689,56661,57100,57587,58300,58926,59543];
        ann.resp.breaths.x2 = [266,826,1311,1929,2499,3137,3596,4307,4539,4945,5151,5494,6048,6608,7172,7714,8268,8827,9460,10152,10637,11327,11881,12524,13010,13553,14234,14735,15203,15831,16544,17087,17562,18185,18800,19296,19840,20383,21122,21586,22220,22846,23468,23991,24255,24814,25379,26007,26416,27055,27699,28163,28801,29435,29989,30530,31084,31639,32187,32736,33443,33995,34465,35172,35700,36280,36908,37373,37672,37920,38310,38600,39339,39825,40442,40922,41469,42102,42662,43210,43833,44308,45013,45546,46174,46718,47082,47409,47974,48675,49159,49792,50199,50816,51454,52014,52489,53025,53669,54371,54856,55722,55970,56268,56685,57092,57588,58358,58912,59546];
    case 15
        ann.resp.breaths.x1 = [620,1154,1533,2067,2607,3054,3485,3852,4383,4883,5152,5791,6261,6817,7265,7746,8207,8559,9072,9450,10011,10446,10993,11517,11930,12648,13509,13948,14461,14913,15415,15885,16380,16850,17389,17976,18402,18822,19365,19778,20300,20757,21248,21865,22335,22676,23293,23750,24767,25333,25911,26430,27009,27513,28057,28570,29113,29617,30128,30585,31146,31637,32080,32637,33124,33633,34157,34743,35252,35691,36091,36591,37004,37550,38076,38533,38950,39480,39972,40520,41076,41591,42109,42630,43148,43670,44152,44639,45215,45641,46237,46680,47267,47615,48917,49213,50400,50830,51348,52811,53254,54224,54841,55350,56857,57400,57743,58296,58778,59526];
        ann.resp.breaths.x2 = [251,599,1148,1522,2050,2604,3016,3491,3879,4391,4882,5141,5832,6243,6803,7251,7693,8173,8542,9080,9350,9909,10474,10980,11516,11933,12662,13015,13511,13954,14440,14968,15425,15905,16364,16865,17361,17868,18385,18821,19370,19797,20330,20758,21259,21792,22325,22698,23284,23759,24207,24751,25273,25811,26437,27023,27540,27984,28585,29065,29651,30087,30583,31132,31617,32098,32678,33079,33628,34196,34776,35257,35679,36064,36586,37014,37550,38125,38553,38927,39450,40020,40521,41075,41559,42118,42582,43115,43659,44224,44593,45166,45652,46211,46712,46939,47272,47673,48079,48932,49217,49813,50341,50906,51333,51966,52389,52740,53310,53811,54234,54455,54867,55373,55727,56228];
    case 16
        ann.resp.breaths.x1 = [407,820,1241,1667,2067,2493,2907,3237,3659,3974,4335,4752,5170,5583,5996,6413,6826,7139,7550,7859,8280,8624,9037,9333,9746,10159,10559,10972,11374,11661,12104,12496,12935,13339,13770,14061,14522,14926,15320,15689,15972,16337,16593,17015,17454,17859,18146,18563,18943,19217,19630,20048,20461,20878,21296,21657,22048,22343,22750,23085,23489,23915,24241,24672,24863,25067,25341,25798,26228,26648,26996,27387,27700,28109,28530,28865,29278,29696,30020,30441,30850,31246,31689,32080,32520,32876,33315,33711,34048,34461,34896,35235,35674,36022,36500,36813,37226,37628,37998,38415,38702,38950,39376,39815,40228,40650,41063,41357,41800,42191,42609,43009,43348,43800,44039,44496,44770,45189,45511,45963,46241,46676,46946,47363,47793,48098,48489,48922,49309,49709,50139,50509,50943,51365,51752,52183,52478,52872,53289,53728,54093,54507,54872,55320,55733,56050,56474,56896,57187,57609,57926,58200,58622,58922,59322,59735];
        ann.resp.breaths.x2 = [409,842,1248,1660,2077,2488,2910,3243,3665,3969,4344,4761,5177,5589,6001,6428,6840,7135,7545,7861,8284,8627,9049,9339,9740,10162,10563,10985,11395,11675,12118,12503,12915,13332,13754,14086,14482,14920,15314,15668,15968,16322,16602,17024,17430,17858,18137,18549,18927,19228,19623,20030,20462,20890,21301,21650,22056,22352,22756,23088,23500,23922,24223,24645,25088,25336,25774,26228,26659,26997,27377,27699,28105,28543,28854,29266,29693,30018,30435,30863,31253,31686,32077,32451,32873,33322,33728,34069,34486,34882,35241,35689,36038,36513,36819,37209,37624,37988,38405,38685,38970,39397,39803,40220,40637,41065,41363,41786,42187,42609,43042,43353,43791,44071,44493,44794];
    case 17
        ann.resp.breaths.x1 = [33,372,663,1050,1467,1867,2259,2693,3163,3541,3943,4287,4683,5100,5543,5974,6343,6809,7274,7689,8015,8359,8707,9276,9867,10337,10650,10941,11270,11665,12087,12630,13135,13570,14022,14348,14852,15176,15593,16041,16454,16841,17311,17641,18076,18433,18809,19183,19583,19922,20409,20830,21217,21652,22009,22417,22889,23489,23933,24350,24750,25467,25998,26500,26970,27448,27843,28274,28630,29165,29626,29930,30280,30715,31198,31650,32163,32576,33085,33472,33896,34248,34735,35070,35548,35974,36491,36974,37461,38007,38411,38872,39372,39828,40289,40767,41211,41800,42283,42783,43243,43700,44226,44687,45159,45628,46211,46593,46976,47407,47898,48354,48628,49074,49357,49757,50135,50583,50987,51496,52043,52576,53076,53524,53963,54493,54933,55441,55859,56378,56817,57326,57839,58513,58930,59309,59639];
        ann.resp.breaths.x2 = [40,356,636,1058,1433,1839,2261,2689,3137,3528,3916,4280,4666,5125,5547,5958,6354,6808,7267,7687,8004,8347,8737,9286,9819,10321,10648,10848,11218,11627,12086,12614,13110,13569,13965,14382,14857,15187,15583,16011,16449,16828,17288,17625,18047,18427,18800,19170,19544,19924,20404,20827,21212,21608,21982,22404,22904,23463,23912,24297,24745,25463,25985,26495,26976,27408,27815,28221,28638,29123,29625,29915,30272,30710,31169,31644,32150,32583,33074,33449,33911,34249,34687,35051,35526,35943,36486,36972,37441,38004,38368,38875,39360,39814,40294,40785,41191,41791,42271,42762,43210,43706,44197,44651,45161,45625,46195,46565,46987,47398,47873,48354,48644,49080,49323,49739];
    case 18
        ann.resp.breaths.x1 = [76,502,911,1328,1750,2167,2598,3007,3428,3835,4257,4670,5087,5496,5917,6326,6743,7196,7580,8011,8407,8837,9237,9654,10076,10515,10924,11322,11761,12143,12587,12978,13387,13835,14265,14657,15093,15507,15907,16346,16728,17176,17593,17998,18415,18839,19257,19661,20078,20509,20904,21339,21752,22174,22593,23015,23433,23833,24233,24680,25093,25515,25937,26365,26748,27178,27583,28017,28409,28878,29257,29691,30107,30515,30915,31363,31780,32185,32620,33024,33437,33839,34283,34696,35113,35522,35930,36343,36752,37183,37602,38028,38433,38867,39272,39689,40107,40511,40941,41348,41761,42217,42604,43022,43439,43891,44283,44691,45089,45528,45911,46385,46780,47176,47633,47980,48446,48843,49304,49696,50104,50513,50948,51396,51800,52183,52593,53015,53428,53893,54237,54685,55120,55502,55963,56335,56804,57200,57617,58061,58435,58861,59283,59696];
        ann.resp.breaths.x2 = [82,499,926,1343,1755,2172,2588,3005,3422,3837,4254,4671,5093,5499,5922,6344,6761,7172,7582,7999,8421,8843,9260,9666,10083,10511,10922,11342,11754,12176,12588,13010,13427,13838,14250,14677,15087,15515,15926,16338,16755,17172,17594,18011,18427,18837,19259,19676,20093,20505,20927,21344,21755,22183,22598,23004,23421,23843,24265,24672,25099,25505,25933,26337,26759,27176,27593,28015,28432,28854,29266,29688,30108,30515,30937,31359,31776,32187,32604,33021,33438,33858,34270,34687,35109,35520,35932,36359,36771,37193,37608,38020,38442,38854,39270,39693,40109,40521,40943,41358,41775,42192,42604,43020,43437,43859,44271,44683,45108,45530,45947,46359,46770,47187,47609,48021];
    case 19
        ann.resp.breaths.x1 = [289,759,1237,1707,2167,2641,3111,3580,3939,4404,4909,5813,6339,7650,8141,8585,9515,9989,10476,10933,11409,11874,12357,12839,13261,13726,14170,14661,15176,15611,16567,17780,18254,18693,19157,19643,20074,20552,21030,21483,21970,22763,23376,23815,24359,24772,25263,25720,26180,26652,27100,27604,28517,29009,29457,29922,30415,30911,31350,31802,32233,32698,34004,34461,34909,35230,35517,35974,36483,36926,37400,37633,38350,38811,39293,39528,40015,40489,40950,41413,41874,42661,43104,43600,43913,44865,45350,45811,46267,46576,47115,47480,48033,48507,48974,49409,49904,50309,50822,51291,51743,52222,52663,53180,53650,54098,54528,54967,55459,55867,56352,56830,57313,57743,58278,58713,59222,59604];
        ann.resp.breaths.x2 = [309,773,1227,1707,2150,2641,3116,3591,3948,4423,4871,5816,6243,6951,7156,7656,8136,8590,9524,10009,10453,10943,11400,11875,12324,12815,13268,13754,14203,14688,15166,15599,16111,16554,16997,17778,18259,18691,19159,19644,20088,20573,21006,21507,21993,22457,22693,23379,23843,24313,24756,25252,25722,26186,26659,27113,27604,27825,28538,28997,29477,28042,29910,30404,30879,31348,31818,32277,32694,33000,33644,33995,34497,34940,35230,35552,36001,36486,36940,37441,37640,38363,38832,39297,39524,40015,40489,40949,41422,41881,42287,42662,43105,43606,43923,44414,44846,45340,45810,46301,47119,47562,48058,48491,48964,49449,49919,50352,50827,51307,51766,52241,52708,53178,53606,54112];
    case 20
        ann.resp.breaths.x1 = [150,580,980,1476,1898,2424,2946,3489,3913,4443,4952,5409,5904,6400,6870,7343,7811,8389,8872,9276,9763,10367,10802,11317,11848,12370,12843,13304,13843,14374,14774,15315,15750,16241,16767,17320,17759,18259,18707,19243,19626,19626,20196,20635,21152,21691,22204,22715,23207,23720,24193,24733,25141,25637,26124,26622,27096,27657,28209,28652,29139,29704,30224,30693,31146,31711,32124,32663,33180,33628,34126,34617,35213,35657,36178,36717,37157,37633,38137,38620,38885,39141,39628,40193,40615,41102,41678,42204,42700,43126,43722,44265,44726,45250,45741,46280,46680,47124,47672,48167,48711,49252,49709,50235,50700,51278,51848,52339,52841,53498,54041,54502,55028,55585,56015,56565,57061,57552,58117,58622,59204,59826];
        ann.resp.breaths.x2 = [145,583,989,1454,1897,2393,2963,3485,3937,4518,4951,5404,5900,6375,6871,7330,7835,8405,8875,9302,9772,10368,10806,11279,11849,12361,12852,13290,13838,14377,14794,15325,15752,16259,16813,17335,17736,18274,18739,19243,19660,20162,20637,21185,21666,22225,22719,23199,23716,24186,24703,25173,25637,26117,26638,27108,27604,28189,28638,29144,29672,30166,30694,31195,31628,32129,32641,33169,33586,34127,34613,35193,35663,36185,36729,37151,37640,38099,38642,39112,39629,40194,40579,41123,41654,42171,42677,43131,43691,44261,44720,45240,45736,46243,46665,47150,47668,48111,48718,49228,49713,50225,50732,51265,51792,52330,52835,53489,54017,54508,55004,55574,56022,56516,57049,57535,58094,58622,59155,59794];
    case 21
        ann.resp.breaths.x1 = [363,885,1415,1963,2480,3046,3454,3835,4187,4665,5043,5500,6013,6704,7274,7654,8063,8450,8880,9320,9672,10046,10467,10854,11250,11674,12165,12635,12991,13500,13874,14330,14796,15198,15611,16050,16537,16950,17407,17954,18398,18809,19230,19604,20017,20487,21235,21696,22043,22676,23037,23467,23937,24372,24815,25346,25798,26220,26726,27126,27604,28391,28852,29204,29678,30046,30489,30833,31163,31507,31928,32476,32950,33398,33943,34257,34643,35004,35400,35809,36209,36578,36991,37365,37767,38207,38576,38993,39380,39793,40202,40537,40941,41283,41691,42039,42383,42739,43126,43630,43952,44374,44774,45128,45528,46007,46407,46863,47228,47659,48024,48472,48800,49217,49613,50043,50396,50813,51243,51539,52004,52391,52785,53150,53480,53846,54198,54576,54911,55328,55785,56124,56452,56878,57165,57478,57813,58139,58509,58883,59187,59570,59913];
        ann.resp.breaths.x2 = [367,884,1427,1971,2493,3032,3454,3832,4196,4666,5046,5505,6006,6687,7299,7661,8062,8468,8891,9328,9661,10051,10479,10864,11258,11675,12150,12651,13020,13495,13870,14345,14804,15193,15620,16058,16544,16945,17420,17953,18417,18827,19233,19623,20024,20499,21249,21703,22077,22650,23041,23463,23959,24355,24809,25315,25774,26228,26706,27102,27609,28042,28406,28838,29213,29677,30050,30504,30858,31190,31491,31939,32478,32963,33401,33943,34275,34613,34993,35404,35827,36201,36565,36977,37378,37766,38173,38563,38985,39381,39814,40178,40537,40917,41268,41706,42039,42366,42730,43147,43611,43949,44371,44788,45140,45536,45995,46427,46865,47230,47657,48021,48443,48816,49212,49608,50056,50383,50832,51243,51523,52014,52367,52788,53131,53479,53864,54191,54582,54925,55331,55790,56133,56464,56881,57171,57472,57815,58147,58506,58870,59197,59567,59920];
    case 22
        ann.resp.breaths.x1 = [67,498,902,1315,1720,2150,2572,2872,3307,3698,4126,4491,4948,5387,5804,6187,6617,6943,7378,7785,8133,8472,8898,9333,9689,10102,10507,10941,11343,11674,12083,12500,12909,13339,13757,14113,14530,14943,15346,15793,16102,16493,16841,17250,17698,18115,18541,18939,19374,19787,20230,20600,21035,21396,21817,22239,22654,23085,23420,23798,24207,24633,25033,25346,25772,26163,26578,26952,27313,27583,28022,28430,28865,29257,29687,30102,30493,30767,31193,31620,32037,32285,32707,33120,33541,33943,34387,34691,35087,35517,35939,36357,36787,37200,37598,38015,38424,38859,39267,39689,39976,40407,40828,41228,41583,41861,42278,42683,43135,43539,43957,44348,44704,45115,45541,45963,46350,46776,47189,47593,48028,48441,48746,49174,49583,50009,50422,50826,51261,51670,52113,52500,52924,53346,53776,54172,54598,55011,55437,55741,56154,56548,57004,57413,57830,58230,58657];
        ann.resp.breaths.x2 = [77,478,894,1327,1739,2156,2578,2873,3290,3702,4117,4491,4940,5367,5800,6191,6613,6951,7362,7793,8120,8463,8901,9323,9672,10094,10511,10933,11353,11675,12092,12498,12925,13327,13765,14108,14535,14941,15367,15768,16106,16512,16855,17261,17689,18100,18522,18927,19349,19761,20193,20610,21027,21402,21813,22230,22656,23062,23421,23801,24207,24613,25025,25331,25759,26165,26590,26944,27303,27593,28020,28432,28833,29255,29667,30098,30504,30773,31195,31617,32024,32282,32704,33116,33538,33937,34359,34692,35104,35520,35932,36344,36776,37188,37587,38009,38431,38859,39270,39693,39988,40394,40811,41234,41564,41870,42287,42693,43110,43527,43944,44345,44693,45103,45536,45953,46359,46770,47187,47609,48021,48438,48744,49175,49592,50009,50431,50827,51265,51681,52093,52513,52930,53336,53769,54175,54598,55004,55442,55737,56165,56564,56986,57398,57830,58205,58659,59007,59414,59841];
    case 23
        ann.resp.breaths.x1 = [224,520,928,1359,1776,2202,2593,3007,3450,3707,4117,4683,5096,5526,5939,6352,6748,7178,7576,8028,8424,8841,9267,9693,10102,10524,10924,11352,11778,12174,12583,13013,13430,13843,14265,14678,15107,15515,15937,16367,16759,17180,17607,18011,18441,18835,19274,19683,20122,20530,20939,21352,21761,22200,22607,23024,23446,23846,24280,24685,25102,25550,25920,26343,26770,27200,27587,28043,28426,28861,29287,29670,30102,30502,31376,31780,32202,32615,32880,33302,33689,34126,34565,34965,35387,35822,36213,36626,37052,37876,38311,38720,39154,39554,39959,40376,40815,41198,41604,42026,42448,42852,43283,43696,44100,44539,44939,45367,45789,46193,46615,47046,47472,47863,48280,48693,49130,49522,49952,50383,50778,51200,51626,52039,52452,52876,53276,53680,54107,54528,54954,55367,55780,56185,56609,57039,57452,57865,58274,58700,59122,59535,59952];
        ann.resp.breaths.x2 = [198,541,926,1348,1770,2193,2604,3011,3433,3702,4154,4703,5088,5526,5922,6359,6761,7188,7587,8025,8426,8838,9260,9693,10109,10526,10933,11337,11775,12171,12598,13010,13443,13849,14266,14677,15103,15520,15937,16359,16760,17193,17609,18026,18443,18842,19280,19687,20135,20520,20948,21365,21761,22193,22592,23020,23447,23859,24286,24672,25115,25521,25938,26363,26749,27192,27577,28031,28437,28865,29276,29672,30124,30504,30963,31354,31786,32182,32620,32879,33290,33707,34122,34539,34966,35378,35795,36222,36629,37046,37872,38305,38706,39133,39540,39978,40368,40811,41207,41611,42028,42445,42857,43279,43706,44086,44535,44947,45367,45784,46195,46617,47024,47441,47852,48280,48691,49159,49518,49945,50362,50784,51185,51613,52019,52452,52856,53284,53685,54133,54518,54962,55352,55801,56186,56622,57018,57456,57862,58284,58691,59113,59535,59952];
    case 24
        ann.resp.breaths.x1 = [98,415,707,1028,1337,1654,1989,2311,2620,2967,3276,3615,3952,4283,4591,4948,5287,5617,5961,6287,6613,6974,7291,7615,7959,8302,8628,8946,9276,9628,9967,10285,10620,10946,11265,11613,11939,12274,12583,12917,13230,13574,13909,14239,14565,14904,15263,15589,15920,16263,16593,16928,17276,17585,17920,18267,18576,18909,19239,19565,19904,20226,20552,20865,21222,21548,21874,22191,22524,22837,23180,23498,23828,24172,24511,24841,25172,25507,25828,26167,26500,26830,27148,27491,27843,28165,28530,28865,29217,29539,29874,30211,30546,30863,31185,31550,31876,32315,32654,32972,33307,33602,33943,34296,34630,34948,35265,35613,35922,36243,36574,36900,37248,37580,37880,38176,38524,38833,39141,39454,39815,40137,40459,40785,41120,41452,41787,42122,42461,42778,43104,43439,43804,44122,44461,44800,45133,45489,45837,46172,46515,46859,47207,47546,47880,48224,48541,48900];
        ann.resp.breaths.x2 = [103,404,715,1026,1348,1644,1976,2303,2620,2963,3280,3612,3958,4286,4613,4940,5267,5615,5953,6286,6608,6966,7309,7629,7951,8294,8621,8938,9265,9629,9962,10289,10595,10954,11274,11606,11939,12261,12598,12915,13242,13580,13907,14245,14561,14904,15251,15588,15921,16253,16591,16918,17266,17594,17921,18253,18575,18900,19233,19565,19903,20225,20547,20863,21217,21534,21866,22193,22518,22840,23168,23489,23822,24170,24508,24846,25173,25500,25827,26165,26511,26817,27155,27498,27846,28168,28527,28859,29203,29535,29878,30214,30546,30863,31195,31544,31908,32303,32646,32984,33301,33623,33948,34301,34634,34945,35267,35605,35927,36249,36571,36903,37235,37582,37877,38168,38521,38822,39128,39471,39782,40136,40447,40790,41128,41458,41775,42118,42466,42778,43105,43437,43796,44134,44451,44794,45140,45483,45821,46169,46522,46871,47214,47557,47889,48222,48538,48900,49249,49592,49935,50283,50605,50943,51286,51650,51982,52315,52698,53109,53458,53785,54139,54466,54803,55115,55437,55764,56096,56427,56738,57081,57419,57746,58100,58427,58759,59086,59424,59757];
    case 25
        ann.resp.breaths.x1 = [63,541,1033,1498,1989,2459,2950,3380,3878,4396,4887,5400,5874,6413,6887,7400,7859,8350,8789,9250,9733,10224,10711,11185,11626,12157,12591,13109,13600,14078,14539,15037,15528,15946,16424,16867,17450,17893,18333,18791,19239,19717,20122,20548,21022,21500,21922,22400,22859,23320,23837,24259,24711,25185,25650,26146,26591,27096,27583,28057,28513,29026,29491,29961,30420,30880,31311,31733,32172,32615,33089,33559,34039,34539,34991,35443,35909,36383,36843,37309,37785,38211,38650,39093,39554,40050,40437,40833,41261,41674,42139,42613,43026,43474,43883,44291,44709,45189,45746,46059,46489,46889,47341,47807,48285,48770,49361,49861,50365,50835,51287,51726,52191,52641,53146,53633,54141,54633,55137,55663,56150,56609,57091,57526,58000,58500,58974,59457,59939];
        ann.resp.breaths.x2 = [71,541,1042,1496,1987,2456,2947,3369,3879,4365,4882,5404,5879,6407,6887,7394,7856,8352,8785,9255,9730,10220,10722,11181,11633,12144,12593,13110,13601,14076,14535,15034,15520,15947,16422,16871,17446,17879,18338,18790,19233,19703,20119,20552,21016,21491,21919,22399,22856,23310,23832,24265,24703,25189,25648,26149,26590,27108,27582,28057,28511,29049,29487,29957,30420,30879,31317,31744,32166,32631,33069,33565,34043,34528,34982,35436,35916,36386,36834,37309,37782,38210,38664,39102,39550,40036,40437,40832,41263,41670,42134,42604,43031,43464,43875,44292,44704,45161,45752,46053,46470,46897,47340,47778,48285,48768,49359,49856,50362,50832,51291,51729,52183,52629,53136,53606,54107,54629,55141,55669,56149,56601,57081,57519,58005,58495,58960,59456,59936];
    case 26
        ann.resp.breaths.x1 = [59,367,1128,1811,2515,2915,3398,3930,4478,4870,5404,5739,6083,6487,6961,7330,7728,7989,8441,8885,9289,9707,10280,10763,11063,11452,11852,12296,12735,13109,13526,13917,14343,14752,15159,15572,15867,16367,16767,17233,17693,18202,18693,19074,19665,20161,20526,20922,21296,21861,22274,22685,23089,23524,23902,24280,24672,25080,25463,25828,26228,26674,27096,27483,27952,28287,28661,29048,29478,29939,30424,30946,31267,31733,32180,32698,33267,33878,34296,34761,35283,35917,36413,36939,37448,38089,38672,39298,39815,40202,40524,41002,41435,42052,42548,43017,43678,44235,44783,45163,45428,45833,46154,46528,46924,47315,47833,48385,48839,49304,50370,51074,51691,52461,52972,53333,53859,54241,54680,55102,55520,55933,56228,56787,57243,57674,58126,58596,59009,59443,59757];
        ann.resp.breaths.x2 = [66,388,583,805,1132,1522,1807,2472,2889,3391,3890,4449,4840,5420,5742,6048,6423,6951,7330,7756,7967,8416,8859,9281,9693,10257,10759,11001,11337,11801,12266,12778,13147,13527,13960,14329,14767,15172,15599,15852,16354,16770,17198,17657,18142,18686,18832,19064,19703,20199,20568,20911,21249,21898,22283,22661,23051,23532,23906,24281,24650,25136,25468,25848,26239,26316,26722,27092,27477,27941,28327,28670,29055,29477,29941,30409,30953,31248,31691,32193,32704,33269,33869,34323,34761,35294,35916,36391,36929,37436,38078,38653,39292,39798,40194,40505,40980,41437,42060,42556,43020,43675,44234,44809,45145,45393,45778,46164,46544,46913,47319,47821,48369,48853,49286,50341,51053,51666,52040,52447,52930,53310,53785,54197,54629,55083,55500,55906,56223,56743,57250,57646,58057,58643,59044,59403,59735];
    case 27
        ann.resp.breaths.x1 = [385,1020,1702,2185,2993,3563,4300,4883,5696,6317,7009,7867,8498,9150,9867,10533,11072,11791,12483,13135,13870,14557,15285,15989,16759,17502,18207,18707,19348,20100,20861,21539,22135,22854,23489,24067,24759,25398,26107,26865,27552,28209,28935,29570,30367,30933,31628,32311,32941,33698,34470,35339,36074,36796,37413,38924,39641,40411,41059,41757,42483,43165,43887,44661,45320,46037,46707,47507,48237,48830,49613,50270,51057,51735,52491,53324,54011,54733,55359,56054,56835,58717];
        ann.resp.breaths.x2 = [388,1058,1697,2182,2942,3644,4317,4940,5626,6381,6977,7846,8500,9112,9851,10447,11144,11823,12498,13184,13944,14625,15282,16016,16828,17509,18206,18879,19381,20177,20863,21608,22209,22883,23484,24149,24840,25516,26112,26954,27630,28263,29013,29614,30383,30963,31675,32409,33016,33863,34650,35478,36175,36877,37561,38421,39028,39714,40474,41263,41875,42625,43163,43918,44667,45409,46111,46802,47478,48153,48827,49666,50357,51053,51750,52619,53389,53985,54751,55220,56102,56643,57556,57915,58290,58701,58997,59245,59461];
    case 28
        ann.resp.breaths.x1 = [241,672,1085,1498,1941,2372,2828,3276,3715,4165,4530,4761,5043,5378,5757,6161,6617,7004,7448,7872,8359,8820,9267,9720,10176,10567,10902,11237,11600,11987,12396,12830,13257,13700,14239,14648,15085,15550,15959,16341,16720,17111,17520,17885,18241,18680,19126,19504,20013,20396,20683,21070,21539,21917,22304,22711,23146,23524,23967,24402,24876,25320,25715,26150,26535,26904,27304,27683,28096,28535,28917,29191,29470,29796,30172,30559,30972,31385,31828,32254,32707,33189,33620,34009,34404,34830,35239,35678];
        ann.resp.breaths.x2 = [235,652,1074,1485,1908,2383,2826,3269,3739,4175,4539,5051,5404,5768,6175,6597,7009,7431,7883,8363,8827,9286,9751,10173,10532,10891,11234,11601,12002,12392,12883,13279,13722,14234,14619,15108,15567,15968,16322,16728,17119,17520,17873,18237,18718,19148,19518,20030,20357,20710,21085,21534,21935,22320,22708,23152,23532,23970,24418,24877,25289,25732,26149,26543,26928,27297,27672,28094,28495,28907,29234,29472,29794,30182,30536,30995,31364,31792,32235,32699,33195,33612,34006,34402,34840,35235,35663,36085,36513,36882,37230,37577,37951,38394,38817,39292,39814,40204,40632,41128,41617,42044,42414,42841,43263,43685,44139,44567,44989,45393,45810,46201,46607,47040,47441,47810,48195,48623,49080,49454,49903,50330,50779,51222,51634,52056,52484,52920,53373,53811,54228,54682,55083,55521,55996,56422,56812,57213,57641,58063,58495,58912,59345,59751];
    case 29
        ann.resp.breaths.x1 = [207,628,1089,1485,1915,2320,2698,3067,3446,3809,4165,4513,4878,5213,5613,6030,6326,6709,7109,7483,7898,8276,8698,9063,9476,9798,10189,10607,11020,11461,11752,12126,12487,12857,13261,13713,14178,14674,15028,15476,15911,16324,16711,17154,17554,17911,18324,18689,19065,19430,19822,20152,20478,20861,21252,21570,21865,22222,22559,22989,23289,23698,24085,24537,24954,25337,25820,26189,26630,27035,27487,27874,28217,28622,29026,29417,29861,30128,30541,30933,31307,31698,32176,32550,32933,33333,33707,34143,34561,34948,35330,35696,36087,36470,36957,37270,37602,38033,38441,38846,39324,39693,40098,40511,40872,41348,41717,42113,42609,43009,43404,43735,44057,44587,45054,45393,45741,46189,46541,47011,47376,47802,48241,48659,49078,49474,49974,50430,50770,51174,51561,51991,52417,52693,53137,53511,53950,54385,54780,55228,55689,56085,56635,57026,57574,58243,58639,59048,59504,59800];
        ann.resp.breaths.x2 = [230,567,1021,1501,1913,2388,2731,3069,3401,3739,4154,4491,4845,5251,5647,5985,6312,6650,7072,7484,7946,8257,8685,9044,9497,9772,10178,10590,11102,11490,11749,12097,12514,12841,13258,13680,14176,14582,15013,15430,15921,16327,16723,17224,17509,17910,18322,18718,19048,19396,19861,20067,20468,20816,21243,21507,21903,22141,22545,22967,23294,23664,24112,24534,24951,25315,25859,26160,26622,27013,27445,27804,28184,28580,28949,29334,29862,30092,30504,30910,31332,31670,32261,32546,32921,33285,33707,34127,34592,34919,35320,35695,36080,36486,36993,37230,37555,38051,38394,38790,39365,39640,40125,40453,40869,41279,41696,42113,42614,42973,43390,43696,44013,44371,44588,45061,45346,45715,46222,46623,46987,47298,47799,48296,48612,49032,49449,49935,50452,50737,51201,51560,52067,52441,52708,53141,53474,53954,54376,54793,55194,55606,56054,56543,57049,57551,58195,58590];
    case 30
        ann.resp.breaths.x1 = [211,680,1133,1498,1898,2350,2854,3259,3672,4109,4483,4861,5248,5643,6017,6396,6722,7087,7572,8093,8641,9102,9507,9959,10359,10754,11128,11487,11922,12352,12778,13109,13513,14057,14596,15054,15515,15915,16341,16754,17141,17511,17893,18254,18659,19017,19435,19796,20191,20526,20874,21252,21657,22161,22672,23089,23476,23854,24241,24641,25050,25563,25946,26339,26752,27213,27622,28043,28461,28787,29217,29661,30033,30420,30893,31320,31728,32076,32515,32885,33272,33641,34065,34474,34870,35252,35635,36026,36413,36800,37200,37598,37928,38289,38702,39085,39498,39907,40280,40685,41037,41378,41778,42191,42583,42878,43239,43657,44030,44435,44813,45193,45589,45985,46350,46746,47354,47772,48254,48620,48983,49396,49817,50217,50609,50978,51335,51713,52083,52448,52837,53233,53563,53950,54337,54724,55093,55463,55837,56198,56535,56922,57252,57657,58030,58365,58617,58930,59417];
        ann.resp.breaths.x2 = [203,678,1121,1522,1892,2351,2852,3269,3681,4111,4486,4892,5257,5647,6022,6396,6713,7088,7571,8104,8642,9107,9513,9946,10358,10764,11128,11501,11912,12350,12788,13094,13506,14065,14604,15050,15525,15931,16348,16749,17145,17530,17889,18264,18660,19022,19439,19792,20172,20531,20863,21259,21676,22141,22682,23099,23479,23854,24249,24624,25046,25595,25964,26358,26759,27213,27625,28052,28474,28791,29229,29656,30018,30420,30894,31322,31723,32071,32509,32889,33248,33649,34069,34491,34861,35257,35631,36027,36428,36797,37199,37603,37925,38299,38701,39102,39492,39904,40289,40685,41059,41395,41791,42192,42582,42873,43237,43648,44034,44424,44815,45187,45604,45974,46369,46734,47351,47789,48248,48617,49006,49396,49829,50235,50605,50974,51338,51703,52093,52452,52856,53215,53574,53959,54339,54740,55094,55495,55838,56191,56548,56923,57282,57672,58047,58369,58596,58965];
    case 31
        ann.resp.breaths.x1 = [124,511,941,1315,1928,2220,2554,3167,3567,3961,4335,4696,5070,5435,5852,6235,6626,7043,7465,7902,8220,8615,8967,9393,9772,10167,10554,10924,11348,11739,12139,12526,12957,13339,13757,14209,14596,14991,15441,15867,16276,16689,17141,17593,17998,18393,19209,19530,19974,20343,20743,21165,21570,22013,22439,22872,23307,23750,24189,24567,24989,25433,25807,26180,26626,26974,27357,27730,28161,28543,28922,29343,29717,30098,30528,30976,31367,31837,32246,32663,33111,33559,33930,34304,34752,35139,35539,35913,36322,36696,37130,37524,37920,38328,38728,39180,39615,39998,40450,40937,41330,41787,42070,42539,43022,43483,43874,44370,44822,45241,45702,46141,46615,47028,47459,47833,48276,48667,49074,49548,49930,50343,50791,51265,51717,52104,52520,52902,53376,53811,54189,54646,55120,55554,55963,56357,56791,57165,57626,58004,58409,58761,59191,59578];
        ann.resp.breaths.x2 = [129,520,931,1322,1913,2562,3185,3570,3953,4323,4697,5056,5452,5842,6233,6644,7046,7478,7898,8210,8616,8975,9397,9782,10141,10558,10917,11358,11733,12144,12519,12957,13327,13754,14203,14572,14978,15441,15858,16264,16697,17135,17599,17995,18391,19191,19534,19961,20357,20747,21164,21571,22014,22457,22856,23310,23743,24191,24566,24983,25421,25806,26186,26627,26986,27356,27725,28163,28532,28928,29340,29699,30113,30551,30958,31385,31839,32235,32668,33121,33575,33937,34312,34755,35156,35552,35927,36333,36729,37098,37524,37914,38352,38727,39181,39608,39983,40437,40949,41342,41780,42055,42535,43031,43495,43870,44366,44830,45235,45673,46132,46596,47040,47483,47847,48274,48670,49064,49523,49956,50346,50795,51259,51718,52104,52513,52898,53379,53811,54197,54656,55109,55569,55943,56369,56770,57181,57582,58026,58400,58775,59208,59577];
    case 32
        ann.resp.breaths.x1 = [85,367,650,933,1211,1533,1833,2111,2389,2741,3011,3307,3602,3861,4161,4457,4791,5065,5370,5696,5987,6291,6613,6904,7213,7524,7828,8115,8441,8746,9072,9350,9676,9959,10254,10580,10907,11198,11478,11791,12078,12374,12639,12970,13304,13630,13909,14196,14474,14761,15041,15389,15672,16002,16298,16615,16889,17198,17520,17820,18111,18420,18720,19048,19374,19670,19970,20287,20565,20870,21139,21435,21700,22052,22335,22624,22963,23250,23524,23815,24124,24446,24763,25085,25363,25680,26002,26296,26583,26843,27143,27422,27713,28043,28339,28622,28909,29261,29517,29813,30128,30454,30759,31098,31376,31746,32028,32311,32589,32902,33220,33550,33883,34165,34474,34813,35148,35452,35739,36017,36352,36665,36943,37226,37515,37789,38120,38398,38702,39046,39337,39685,39993,40263,40554,40850,41111,41391,41683,42026,42361,42661,42952,43270,43561,43896,44191,44513,44813,45093,45402,45693,45989,46302,46628,46937,47220,47563,47872,48172,48446,48770,49104,49400,49704,49996,50348,50665,50965,51322,51622,51930,52226,52541,52837,53146,53498,53807,54098,54428,54737,55028,55346,55680,55989,56300,56561,56870,57183,57496,57791,58117,58443,58770,59074,59357,59722];
        ann.resp.breaths.x2 = [82,372,646,926,1201,1538,1828,2119,2393,2736,3016,3317,3596,3869,4148,4433,4787,5077,5352,5689,5995,6275,6618,6914,7193,7529,7825,8109,8453,8753,9038,9371,9677,9962,10247,10584,10891,11170,11501,11780,12086,12361,12641,12973,13295,13627,13902,14208,14466,14751,15040,15377,15673,15974,16306,16612,16887,17166,17520,17810,18095,18427,18728,19016,19375,19650,19987,20283,20568,20853,21127,21407,21687,22051,22330,22608,22962,23257,23526,23811,24139,24450,24745,25088,25363,25669,26001,26284,26580,26854,27134,27424,27699,28057,28337,28617,28896,29250,29530,29809,30092,30451,30741,31090,31380,31712,32024,32309,32588,32868,33232,33522,33885,34164,34460,34808,35151,35452,35742,36016,36365,36660,36935,37220,37508,37777,38120,38410,38690,39044,39339,39687,39988,40268,40547,40832,41102,41379,41670,42013,42366,42651,42947,43242,43522,43881,44166,44514];
    case 33
        ann.resp.breaths.x1 = [102,620,1180,1689,2198,2720,3293,3711,3826,4435,5039,5600,6135,6739,7339,8093,8663,9176,9772,10320,10815,11400,11957,12504,13109,13604,14078,14630,15289,15867,16593,17159,17902,18450,19057,19570,20174,20761,21409,22000,22589,23193,23598,24189,24672,25202,25802,26391,26978,27404,28226,28961,29830,30511,30954,31359,31807,32307,32824,33424,34026,34557,35117,35678,36265,36783,37261,37763,38263,38863,39433,39959,40554,41180,41661,42157,42670,43139,43635,44165,44726,45220,45715,46176,46650,47176,47811,48341,48904,49378,49970,50465,51191,51743,52400,52872,53380,53859,54380,54898,55415,55959,56565,57230,57983,58652,59326];
        ann.resp.breaths.x2 = [98,631,1179,1691,2198,2720,3296,3837,4439,5035,5610,6154,6782,7462,8088,8658,9207,9772,10321,10822,11395,11949,12493,13089,13633,14097,14672,15288,15905,16617,17235,17900,18485,19059,19586,20193,20790,21465,22046,22603,23231,23616,24202,24693,25247,25848,26406,27002,27561,28242,28976,29857,30499,31000,31401,31834,32346,32858,33459,34048,34618,35172,35716,36270,36782,37251,37756,38273,38885,39439,40009,40584,41202,41643,42150,42662,43126,43659,44155,44709,45219,45704,46174,46686,47256,47784,48338,48911,49439,49972,50505,51233,51713,52389,53373,53843,54376,54914,55416,56012,56617,57319,57962,58643,59319,59994];
    case 34
        ann.resp.breaths.x1 = [293,759,1089,1472,2102,2724,2950,3363,3737,4130,4517,4913,5352,5765,6239,6696,7152,7559,7828,8224,8733,9063,9559,9876,10346,10780,11167,11861,12291,13926,14352,14678,15093,15372,15698,16046,16385,16702,17041,17511,17928,18376,18843,19252,19717,20078,20487,20965,21383,21739,22161,22528,22893,23280,23676,24033,24450,24941,25398,25820,26270,26691,27122,27583,27939,28387,28809,29239,29674,30159,30667,31128,31572,32024,32454,32902,33337,33728,34083,34443,34848,35291,35700,36161,36683,37187,37650,38085,38520,38954,39428,39863,40337,40776,41211,41665,42122,42548,43017,43461,43922,44391,44817,45211,45611,46033,46428,46863,47263,47746,48180,48715,49187,49626,50091,50513,50913,51365,51848,52348,52859,53380,53872,54341,54776,55220,55693,56133,56583,57026,57435,57870,58330,58839,59304,59765];
        ann.resp.breaths.x2 = [298,763,1106,1485,2087,2715,2958,3348,3734,4133,4523,4914,5346,5763,6238,6697,7156,7566,7840,8210,8716,9054,9540,9877,10331,10727,11149,11327,11865,12282,12593,12952,13173,13385,13923,14366,14688,15119,15383,15694,16053,16391,16718,17055,17525,17942,18396,18848,19254,19734,20088,20484,20958,21391,21766,22151,22529,22893,23268,23669,24049,24466,24956,25400,25822,26263,26696,27118,27577,27947,28385,28823,29245,29677,30161,30668,31121,31586,32024,32472,32921,33332,33734,34101,34454,34845,35294,35726,36180,36681,37183,37661,38094,38511,38980,39418,39872,40342,40769,41207,41670,42123,42551,43026,43453,43918,44382,44815,45208,45620,46026,46449];
    case 35
        ann.resp.breaths.x1 = [350,820,1172,1620,1920,2298,2641,3107,3589,4000,4435,4617,4835,5109,5417,5783,6283,6687,7000,7478,7924,8393,8893,9354,9798,10254,10693,11133,11604,12091,12600,13000,13504,14004,14439,14904,15341,15802,16263,16667,17085,17546,17980,18411,18830,19291,19726,20165,20483,20883,21248,21765,22257,22772,23276,23733,24185,24415,24893,25337,25789,26267,26713,27152,27609,28035,28422,28817,29222,29687,30124,30498,30937,31307,31763,32193,32580,32959,33328,33680,34157,34596,35065,35496,35904,36348,36783,37200,37693,38124,38620,39028,39385,39715,40120,40546,41015,41452,41917,42343,42726,43165,43548,44000,44439,44822,45254,45693,46159,46585,47080,47528,47985,48415,48741,49217,49639,49952,50122,50561,51009,51404,51787,52243,52667,53111,53533,53954,54328,54707,55033,55437,55889,56330,56796,57152,57513,57943,58378,58726];
        ann.resp.breaths.x2 = [377,794,1164,1612,1913,2298,2631,3100,3565,4016,4439,4856,5093,5431,5763,6265,6687,7003,7484,7920,8384,8875,9344,9809,10263,10685,11133,11590,12086,12572,13042,13495,13981,14424,14883,15335,15810,16269,16675,17103,17525,17995,18417,18832,19270,19713,20162,20499,20863,21265,21766,22267,22756,23252,23701,24165,24408,24872,25331,25785,26234,26696,27155,27609,28036,28416,28791,29255,29683,30119,30493,30905,31332,31776,32193,32583,32942,33338,33681,34164,34576,35040,35478,35900,36333,36766,37209,37687,38152,38611,39049,39376,39714,40109,40558,41007,41443,41912,42340,42746,43184,43553,43991,44424,44820,45245,45694,46158,46596,47066,47525,47984,48427,48790,49228,49660,50141,50573,51006,51381,51808,52225,52666,53109,53558,53943,54307,54714,55025,55431,55891,56316,56791,57176,57503,57941,58379,58759,59171,59704];
    case 36
        ann.resp.breaths.x1 = [454,867,1259,1646,2067,2446,2820,3224,3672,4091,4483,4878,5313,5717,6043,6230,6552,6948,7387,7767,8163,8576,9007,9398,9837,10250,10659,11093,11496,11896,12291,12570,12943,13413,13678,14052,14435,14826,15202,15680,16089,16493,16859,17263,17693,18063,18428,18852,19252,19696,20083,20435,20830,21148,21548,21974,22374,22741,23141,23593,23967,24372,24807,25154,25515,25946,26361,26726,27135,27552,27922,28357,28717,29130,29500,29926,30350,30754,31107,31511,31985,32393,32798,33211,33576,33978,34322,34639,35070,35452,35883,36300,36730,37183,37637,38093,38528,38933,39385,39815,40215,40615,41028,41483,41861,42257,42722,43187,43543,43978,44383,44739,45137,45572,46154,46554,47020,47354,47746,48211,48585,49043,49465,49883,50265,50709,51113,51522,51848,52265,52680,53102,53546,54015,54354,54650,55111,55633,56033,56465,56843,57274,57665,58087,58539,58957,59348,59765];
        ann.resp.breaths.x2 = [478,894,1227,1617,2103,2462,2826,3211,3702,4122,4528,4929,5336,5747,6048,6581,6977,7394,7798,8199,8616,9033,9418,9846,10257,10690,11128,11490,11886,12287,12546,12957,13427,13654,14060,14429,14841,15256,15683,16121,16501,16855,17277,17683,18095,18438,18869,19286,19692,20098,20431,20827,21154,21544,21977,22378,22751,23157,23606,23996,24392,24819,25136,25511,26001,26406,26770,27155,27572,27973,28332,28717,29129,29535,29952,30377,30741,31100,31491,31987,32409,32815,33222,33586,33985,34323,34639,35061,35457,35911,36333,36755,37188,37635,38104,38537,38938,39371,39803,40204,40611,41022,41485,41838,42245,42725,43147,43543,43997,44361,44741,45156,45573,46169,46565,47018,47372,47768,48216,48570,49053,49465,49877,50294,50705,51117,51534,51850,52272,52687,53136,53569,54044,54365,54693,55115,55637,56054,56464,56875,57276,57688,58110,58538,58960];
    case 37
        ann.resp.breaths.x1 = [59,502,924,1324,1846,2280,2702,3080,3533,3948,4270,4800,5313,5983,6435,6887,7317,7741,8093,8446,8789,9133,9480,9911,10228,10641,11085,11487,11817,12248,12600,12987,13457,13796,14287,14543,14970,15402,15811,16298,16724,17085,17511,17941,18167,18720,19126,19726,20096,20422,20843,21339,21709,22091,22550,22937,23320,23759,24141,24498,24876,25280,25654,26050,26457,26874,27217,27657,28074,28530,28926,29357,29713,30093,30407,30802,31146,31524,31933,32215,32641,33002,33393,33737,34122,34513,34939,35291,35635,36009,36409,36683,37165,37567,37902,38389,38798,39185,40028,40454,40963,41317,41683,42083,42417,42809,43204,43578,43926,44270,44709,45146,45476,45759,46289,46685,47089,47367,47776,48176,48498,48852,49248,49691,49983,50365,50839,51374,51796,52191,52602,52985,53333,53585,53776,54154,54737,55267,55693,56167,56570,56965,57335,57774,58135,58465,58974,59400,59830];
        ann.resp.breaths.x2 = [61,520,931,1343,1850,2272,2689,3074,3522,3958,4296,4782,5299,5974,6454,6908,7325,7735,8083,8426,8780,9144,9476,9888,10199,10637,11080,11474,11823,12245,12593,12984,13437,13786,14276,14546,14973,15404,15821,16290,16728,17077,17493,17910,18169,18712,19122,19713,20114,20447,20879,21328,21718,22093,22545,22930,23289,23769,24112,24476,24856,25257,25653,26044,26406,26907,27213,27635,28068,28485,28907,29319,29662,30077,30409,30741,31158,31559,31897,32230,32625,32995,33396,33805,34138,34470,34892,35309,35637,36053,36386,36660,37141,37566,37904,38379,38769,39149,40046,40474,40906,41327,41717,42081,42419,42772,43237,43596,43933,44271,44693,45135,45446,45757,46280,46712,47034,47372,47794,48185,48459,48800,49212,49708,49977,50362,50863,51370,51782,52188,52598,53030,53331,53595,53774,54154,54766,55252,55690,56154,56548,56954,57319,57757,58137,58458,58970,59392,59820];
    case 38
        ann.resp.breaths.x1 = [302,763,1307,1889,2407,2959,3459,4004,4487,4978,5496,5978,6426,6943,7413,7859,8250,8780,9267,9733,10254,10698,11120,11583,12070,12561,13004,13504,13952,14374,14804,15328,15780,16233,16759,17211,17711,18180,18663,19091,19596,20057,20483,20965,21452,22022,22439,23480,24037,24528,24889,25780,26085,26448,26817,27187,27535,27804,28174,28565,29004,29309,29757,30150,30520,30763,31089,31641,31976,32485,32976,33459,33887,34378,34874,35296,35839,36304,36817,37235,37737,38133,38602,39133,39563,40120,40633,41124,41609,42074,42409,42800,43330,43865,44343,44809,45354,45867,46459,46967,47446,48011,48502,49083,50135,50548,50930,51261,51713,52126,52976,53550,54537,55137,56041,56652,57070,57661,58087,58691,59178,59761];
        ann.resp.breaths.x2 = [309,778,1306,1897,2420,2958,3470,4016,4497,4998,5515,5985,6449,6935,7431,7877,8252,8790,9276,9740,10226,10674,11107,11580,12071,12556,12994,13490,13907,14350,14794,15346,15778,16232,16755,17214,17699,18179,18633,19085,19565,20035,20499,20985,21444,22003,22457,22640,22967,23495,24017,24455,24872,25442,25780,26091,26448,26823,27203,27540,27788,28173,28564,28949,29287,29772,30156,30509,31106,31586,32003,32493,32974,33475,33911,34402,34850,35341,35842,36286,36776,37220,37703,38120,38600,39117,39550,40125,40653,41128,41596,42086,42408,42836,43258,43875,44261,44788,45361,45873,46475,47003,47435,47995,48522,49212,50125,50505,50932,51265,51718,52104,52555,52993,53574,54049,54228,54566,55104,55431,56102,56606,57044,57630,58063,58643,59171,59772];
    case 39
        ann.resp.breaths.x1 = [128,650,1141,1728,2172,2685,3241,3693,4191,4704,5204,5817,6313,6878,7404,7967,8450,8937,9389,9876,10315,10946,11409,11861,12326,12939,13478,14026,14552,15050,15572,16046,16533,17015,17563,18111,18646,19165,19717,20204,20717,21261,21757,22248,22772,23302,23785,24272,24815,25363,25859,26343,26887,27422,27943,28452,28970,29452,29948,30433,30946,31476,32007,32537,33046,33611,34091,34613,35122,35661,36187,36735,37270,37793,38307,38811,39363,39837,40346,40811,41224,41757,42283,42835,43235,43730,44209,44617,45120,45585,46015,46567,47137,47746,48393,49009,49774,50426,51043,51713,52265,52841,53415,53989,54533,55128,55772,56330,56843,57461,57996,58487,59057,59578];
        ann.resp.breaths.x2 = [140,657,1153,1723,2193,2689,3211,3686,4201,4708,5241,5811,6323,6845,7425,7972,8447,8896,9381,9851,10336,10922,11400,11875,12361,12931,13474,14034,14535,15055,15573,16079,16554,17040,17583,18121,18644,19164,19724,20204,20737,21249,21734,22246,22788,23315,23806,24313,24830,25347,25859,26358,26907,27435,27957,28458,28970,29445,29962,30430,30958,31485,32018,32525,33058,33612,34090,34618,35141,35658,36191,36734,37262,37766,38294,38806,39355,39856,40352,40817,41228,41775,42271,42825,43247,43728,44224,44604,45098,45583,46037,46565,47135,47726,48375,49016,49766,50436,51059,51697,52262,52835,53400,53985,54550,55136,55774,56332,56854,57461,57984,58485,59028,59577];
    case 40
        ann.resp.breaths.x1 = [89,602,1085,1576,2115,2602,3046,3385,3724,4148,4609,5022,5574,6070,6552,7022,7533,8093,8602,9098,9485,9872,10328,10780,11189,12617,13170,13670,14174,14630,15107,15589,16059,16480,16967,17393,17920,18354,18559,18804,19087,19422,19574,20083,20722,21178,21422,22004,22487,22924,23337,23802,24267,24720,25211,25520,26080,26552,26996,27461,27870,28335,28578,29274,29665,30667,31120,31580,32128,32698,33159,33913,34417,34961,35470,35891,36848,37754,38180,38680,39289,39850,40359,40854,41500,41961,42461,42874,43357,43861,44500,44948,45407,45854,46480,46733,46976,47363,47780,48207,48715,49348,49813,50252,50609,51074,51513,51730,51970,52487,53011,53467,53972,54446,55020,55550,56015,56500,56909,57309,57626,57852,58422,58900,59435,59952];
        ann.resp.breaths.x2 = [66,615,1079,1554,2108,2604,3016,3406,3739,4159,4629,5019,5542,6043,6560,7035,7540,8104,8600,9086,9492,9883,10342,10764,11175,11453,11722,11954,12619,13158,13691,14160,14635,15108,15588,16047,16485,16945,17377,17931,18369,18832,19080,19433,20077,20710,21191,22030,22484,22935,23336,23796,24270,24703,25183,25547,26107,26532,27002,27451,27867,28279,28596,29023,29635,30620,31127,31602,32129,32694,33121,33369,33932,34423,34956,35441,35911,36412,36850,37246,37724,38162,38658,38949,39281,39867,40368,40880,41175,41527,41981,42472,42878,43363,43859,44509,44952,45420,45847,46280,46512,46728,46982,47367,47805,48201,48707,49338,49782,50235,50600,51085,51513,51729,51987,52489,53009,53600,53975,54455];
    case 41
        ann.resp.breaths.x1 = [141,598,1007,1550,2254,2772,3385,3930,4470,4939,5474,5987,6570,7083,7715,8367,9080,9680,10389,10876,11470,11978,12478,13122,13657,14117,14657,15159,15711,16420,16954,17380,17941,18585,19100,19600,20070,20670,21230,21683,22252,22750,23367,23924,24433,24907,25515,26020,26635,27113,27674,28187,28717,29287,29861,30315,30846,31454,31998,32459,32976,33485,34030,34526,34891,35496,36235,36726,37209,37820,38285,38746,39272,39785,40346,41220,41757,42352,42870,43361,43970,44439,44883,45333,45798,46376,46928,47511,48020,48637,49204,49674,50139,50604,51161,51543,52009,53072,53585,54198,54754,55376,56050,56626,57239,58070,58670,59096,59683];
        ann.resp.breaths.x2 = [135,583,1005,1565,2245,2794,3375,3953,4460,4935,5499,5990,6555,7114,7714,8352,9065,9656,10410,10896,11527,12023,12519,13131,13659,14144,14651,15182,15694,16475,16960,17404,17942,18607,19111,19602,20093,20658,21233,21692,22230,22777,23363,23917,24455,24920,25516,26044,26638,27139,27683,28189,28722,29303,29862,30330,30831,31480,31997,32462,32974,33470,34027,34528,35536,36206,36729,37199,37809,38299,38753,39260,39788,40368,41228,41749,42340,42862,43395,43949,44424,44889,45340,45831,46385,46950,47530,48021,48660,49196,49692,50130,50600,51154,51528,52019,52513,53073,53600,54218,54761,55368,56049,56627,57224,58078,58670,59049,59704];
    case 42
        ann.resp.breaths.x1 = [185,624,1076,1541,1967,2463,2872,3293,3746,4196,4639,5048,5465,5891,6339,6787,7204,7641,8063,8511,8915,9354,9772,10189,10628,11041,11461,11900,12309,12765,13170,13583,14017,14422,14857,15341,15785,16250,16702,17098,17559,18046,18498,18991,19417,19870,20270,20670,21109,21522,21943,22326,22746,23167,23641,24050,24480,24920,25324,25733,26185,26604,27022,27439,27922,28400,28796,29283,29687,30141,30576,31002,31433,31837,32307,32737,33141,33585,34061,34517,34952,35409,35848,36304,36778,37174,37637,38080,38476,38889,39333,39789,40237,40667,41059,41487,41922,42387,42839,43265,43748,44157,44561,44978,45385,45785,46207,46633,47037,47472,47920,48324,48724,49148,49574,50022,50474,50913,51387,51861,52261,52707,53133,53541,53976,54441,54937,55380,55780,56172,56583,56983,57400,57822,58291,58748,59200,59643];
        ann.resp.breaths.x2 = [172,620,1090,1538,1976,2478,2873,3317,3734,4217,4660,5056,5468,5895,6328,6776,7214,7645,8062,8516,8922,9350,9782,10189,10621,11044,11474,11902,12303,12762,13179,13585,14013,14445,14867,15340,15789,16253,16712,17103,17562,18063,18496,18974,19418,19882,20278,20695,21111,21518,21961,22330,22745,23168,23637,24065,24492,24914,25321,25748,26186,26601,27018,27461,27920,28395,28817,29271,29709,30140,30578,31021,31443,31844,32303,32747,33164,33591,34064,34518,34961,35415,35837,36312,36776,37172,37635,38067,38479,38922,39344,39782,40241,40669,41075,41480,41944,42403,42857,43274,43717,44155,44577,44973,45383,45773,46216,46644,47045,47483,47926,48322,48734,49133,49571,50024];
    case 43
        ann.resp.breaths.x1 = [107,567,985,1450,1924,2433,2867,3320,3830,4252,4683,5100,5552,5926,6357,6778,7235,7611,8007,8776,9233,9672,10124,10502,10907,11287,11778,12257,12674,13204,13622,14070,14552,14983,15415,15798,16233,16667,17063,17498,17928,18393,18822,19239,19730,20213,20626,21109,21522,21961,22417,22859,23285,23672,24098,24528,24967,25367,25811,26233,26674,27096,27513,27983,28443,28861,29261,29691,30111,30502,30941,31433,31915,32372,32850,33280,33685,34143,34600,35039,35478,35935,36361,36817,37261,37667,38176,38602,39054,39489,39933,40376,40776,41211,41635,42052,42465,42891,43291,43700,44157,44561,44952,45402,45833,46263,46676,47041,47407,47759,48172,48537,48891,49291,49648,50022,50448,50852,51261,51670,52070,52452,52902,53307,53720,54098,54507,54928,55346,55750,56141,56565,57017,57417,57857,58278,58704,59122,59561,59961];
        ann.resp.breaths.x2 = [124,562,989,1464,1934,2430,2879,3322,3821,4265,4692,5109,5552,5948,6359,6776,7257,7619,8009,8790,9228,9666,10115,10532,10901,11305,11770,12266,12667,13195,13622,14060,14556,14978,15420,15826,16222,16654,17050,17488,17968,18401,18821,19228,19745,20220,20615,21122,21534,21993,22410,22888,23299,23695,24102,24529,24946,25358,25817,26279,26696,27118,27530,27989,28453,28838,29261,29688,30113,30515,30953,31433,31923,32388,32842,33269,33718,34148,34602,35061,35457,35948,36386,36813,37262,37672,38204,38600,39075,39503,39930,40363,40790,41239,41670,42060,42477,42889,43311,43696,44160,44556,44968,45388,45842,46264,46702,47045,47409,47763,48195,48549,48900,49301,49660,50035,50436,50837,51265,51655,52082,52457,52893,53315,53722,54107,54524,54930,55363,55759,56139,56569,57028,57435,57862,58311,58701,59118,59561,59968];
    case 44
        ann.resp.breaths.x1 = [393,802,1250,1720,2185,2624,3037,3463,3935,4326,4774,5217,5678,6148,6565,7030,7439,7872,8328,8785,9233,9663,10133,10554,10972,11443,11826,12252,12774,13165,13609,14113,14548,15015,15428,15893,16285,16724,17150,17611,18063,18541,18991,19461,19930,20352,20813,21243,21691,22161,22689,23028,23480,23924,24389,24828,25302,25724,26185,26648,27100,27526,28048,28474,28970,29448,29896,30346,30824,31311,31728,32211,32659,33141,33576,33991,34457,34917,35330,35778,36187,36678,37087,37580,38020,38476,38915,39376,39824,40276,40728,41172,41657,42157,42591,43043,43483,43961,44404,44865,45280,45759,46228,46680,47207,47624,48063,48511,48983,49417,49887,50370,50843,51304,51787,52183,52654,53141,53589,54033,54515,55015,55493,55941,56426,56887,57374,58039,58439,58900,59339,59757];
        ann.resp.breaths.x2 = [298,625,1037,1617,1923,2430,2894,3253,3723,4154,4544,5135,5468,6038,6433,6829,7341,7672,8088,8616,9065,9555,10020,10373,10790,11202,11685,12081,12503,12984,13385,13875,14297,14894,15372,15752,16190,16522,17040,17430,17889,18327,18842,19317,19824,20257,20652,21138,21534,22014,22478,22935,23358,23796,24207,24708,25199,25574,26033,26516,26944,27440,27852,28369,28775,29340,29767,30193,30694,31116,31617,32040,32536,33032,33454,33879,34354,34792,35204,35679,36106,36539,37003,37468,37925,38352,38832,39244,39708,40083,40637,41054,41474,42034,42487,42889,43337,43849,44329,44767,45177,45625,46137,46575,47050,47478,47926,48412,48863,49301,49708,50241,50742,51170,51687,52030,52513,52978,53489,53906,54392,54846,55368,55843,56300,56765,57245,57641,58332,58738,59224,59651];
    case 45
        ann.resp.breaths.x1 = [367,780,1263,1746,2280,2754,3207,3698,4209,4657,5148,5626,6091,6552,7087,7572,8007,8493,9041,9541,10080,10507,10980,11443,11913,12400,12878,13300,13787,14278,14700,15185,15646,16337,16767,17220,17733,18189,18667,19130,19622,20048,20530,21004,21465,21909,22378,22802,23233,23680,24167,24607,25085,25593,26033,26483,26974,27413,27952,28391,28883,29330,29800,30285,30763,31254,31715,32211,32667,33141,33602,34096,34552,34983,35435,35874,36317,36796,37270,37667,38141,38585,39024,39480,39885,40341,40785,41215,41670,42143,42613,43078,43583,44052,44483,44930,45398,45802,46272,46767,47180,47654,48133,48607,49043,49470,49991,50430,50904,51326,51809,52313,52793,53250,53741,54207,54680,55150,55598,56076,56504,57000,57439,57861,58309,58787,59278,59739];
        ann.resp.breaths.x2 = [409,757,1227,1744,2303,2757,3243,3758,4228,4634,5151,5610,6101,6549,7082,7571,8004,8505,9044,9545,10083,10558,11038,11422,11960,12387,12830,13279,13812,14276,14704,15219,15668,16359,16760,17198,17789,18232,18665,19101,19660,20093,20578,21006,21418,21956,22373,22814,23263,23637,24181,24561,25125,25584,26038,26485,27002,27419,27947,28358,28875,29303,29830,30288,30799,31243,31749,32256,32720,33179,33607,34064,34555,34993,35431,35890,36344,36824,37272,37708,38173,38584,39017,39450,39867,40358,40774,41268,41691,42144,42635,43078,43575,43997,44509,44941,45383,45794,46301,46792,47208,47683,48148,48617,49048,49465,49993,50389,50900,51323,51829,52336,52788,53247,53732,54197,54666,55157,55606,56049,56501,56997,57414,57846,58290,58828,59297,59751];
    case 46
        ann.resp.breaths.x1 = [4422,4965,6000,6752,8489,8933,9385,10733,11059,12917,13378,13852,14343,14939,15785,16237,16693,17380,17798,18359,18870,19452,20430,20857,22243,22902,25111,26970,27600,29204,29483,29909,30398,30724,31093,31589,32224,32637,33054,33459,33922,34483,34926,35457,35978,36413,36874,37074,37452,37589,38015,38467,38946,39502,39828,40193,40780,41089,41370,42309,42891,43870,44343,44839,45333,45807,46207,46667,47215,47654,48211,48689,49187,49613,50004,50535,51057,51522,52017,52465,53076,53554,54067,54524,54985,55480,55924,56361,56813,57309,57896,58409,58917,59413,59891];
        ann.resp.breaths.x2 = [309,1005,1602,2113,2377,2678,3190,3580,3734,4454,5003,5673,5985,6243,6734,7199,7730,8489,8954,9397,9978,10273,10722,11022,11606,12424,12857,13427,13902,14392,14978,15810,16216,16691,17361,17815,18391,18853,19460,20473,20837,21444,21892,22241,22877,23547,24197,24814,25104,25352,25669,25869,26086,26564,26976,27530,28432,29187,29487,29947,30319,30752,31174,31665,32156,32673,33063,33464,33906,34391,34982,35436,36006,36449,36877,37077,37494,38057,38484,38970,39545,39846,40183,40600,40817,41218,41902,42287,42947,43875,44382,44830,45361,45794,46211,46686,47208,47694,48206,48712,49185,49613,50030,50557,50985,51528,52030,52473,53062,53584,54059,54524,55009,55505,55949,56374,56838,57324,57915,58416,58923,59424,59904];
    case 47
        ann.resp.breaths.x1 = [133,476,824,1141,1498,1846,2198,2520,2837,3193,3546,3883,4209,4557,4904,5291,5613,5796,5957,6217,6443,6770,7122,7461,7811,8150,8489,8828,9180,9520,9841,10185,10546,10885,11185,11535,12065,12361,12578,12874,13161,13500,13835,14178,14530,14865,15207,15546,15885,16228,16572,16920,17254,17585,17933,18276,18620,18952,19287,19635,19978,20326,20661,20987,21348,21687,22026,22357,22698,23046,23389,23750,24067,24407,24750,25102,25441,25772,26115,26452,26796,27135,27470,27813,28152,28500,28848,29178,29526,29865,30207,30537,30889,31228,31567,31907,32250,32593,32941,33280];
        ann.resp.breaths.x2 = [135,483,815,1142,1491,1844,2182,2530,2847,3195,3544,3885,4201,4560,4903,5283,5600,5958,6265,6449,6766,7119,7468,7809,8146,8484,8832,9186,9518,9872,10189,10542,10885,11186,11538,12060,12577,12867,13168,13490,13838,14181,14524,14867,15214,15551,15889,16232,16570,16908,17251,17578,17931,18274,18623,18953,19291,19639,19987,20315,20663,20990,21338,21692,22024,22357,22703,23057,23389,23732,24059,24408,24751,25088,25442,25769,26107,26458,26791,27144,27466,27809,28158,28506,28838,29176,29519,29862,30214,30536,30884,31232,31575,31913,32240,32588,32937,33274,33612,33948,34291,34634,34987,35309,35658,35995,36349,36687,37014,37362,37708,38051,38384,38716,39065,39413,39751,40104,40421,40769,41112,41458,41796,42123,42472,42820,43158,43501,43828,44187,44524,44867,45203,45536,45879,46232,46575,46908,47240,47594,47931,48280,48602,48953,49301,49644,49966,50315,50658,51006,51338,51676,52024,52367,52703,53051,53389,53732,54054,54408,54756,55094,55437,55769,56112,56464,56801,57123,57472,57820,58163,58506,58833,59176,59530,59867];
    case 48
        ann.resp.breaths.x1 = [520,867,1185,1537,1993,2528,3020,4530,5061,5852,6370,6822,7291,7741,8128,8563,9028,9507,10011,10472,10959,11387,11835,12278,12726,13196,13696,14165,14604,14996,15411,15833,16293,16702,17137,17554,17976,18411,18787,19161,19565,20009,20430,20865,21291,21748,22143,22515,22915,23354,23763,24163,24628,25028,25472,25924,26330,26739,27191,27613,28043,28500,28917,29357,29917,30398,30828,31254,31663,32141,32546,32976,33376,33757,34230,34596,34909,35257,35626,36057,36496,36913,37383,37780,38202,38637,39050,39502,39902,40285,40641,41015,41452,41852,42248,42635,43030,43422,43826,44226,44613,45037,45441,45815,46224,46602,46972,47350,47715,48098,48502,48883,49278,49674,50039,50370,50717,51078,51504,51883,52265,52659,53002,53398,53798,54146,54615,55328,55872,55993,56172,56870,57313,57683,58083,58452,58752,59061,59413,59722];
        ann.resp.breaths.x2 = [77,478,863,1185,1533,1992,2541,3032,3596,3900,4544,5046,5885,6312,6808,7294,7735,8120,8569,9033,9508,10009,10468,10964,11379,11849,12276,12725,13195,13712,14160,14609,14994,15425,15821,16280,16702,17135,17546,17984,18417,18795,19180,19581,20024,20425,20874,21301,21745,22130,22518,22893,23363,23748,24186,24629,25051,25468,25917,26316,26759,27197,27635,28052,28516,28907,29371,29910,30377,30836,31237,31675,32145,32530,32974,33364,33763,34238,34597,34914,35251,35605,36064,36502,36903,37378,37782,38189,38632,39049,39492,39904,40299,40658,41007,41437,41838,42250,42635,43020,43416,43838,44239,44609,45034,45425,45826,46206,46612,46982,47356,47726,48095,48491,48885,49296,49666,50035,50362,50721,51085,51491,51871,52283,52650,52999,53384,53774,54186,54613,55162,55991,56875,57319,57709,58084,58432,58733,59060,59414,59704];
    case 49
        ann.resp.breaths.x1 = [537,841,1302,1754,2163,2641,3076,3498,4009,4443,4870,5287,5891,6296,6770,7235,7711,8111,8576,9093,9567,10033,10511,10980,11404,11896,12387,12874,13326,13817,14304,14743,15237,15685,16111,16524,16963,17354,17767,18159,18602,19091,19535,20000,20296,20526,21030,21239,21548,22078,22524,22937,23333,23715,24137,24520,24893,25150,25480,25911,26309,26709,27148,27561,27983,28157,28483,28917,29343,29822,30176,30685,31215,31646,32063,32472,32924,33415,33939,34461,35000,35548,35961,36400,36948,37326,37724,38137,38502,38759,39133,39550,39954,40385,40802,41233,41683,42122,42578,43000,43417,43874,44252,44696,45107,45528,45924,46420,46837,47259,47767,48241,48676,49143,49509,49904,50322,50813,51239,51757,52187,52659,53093,53467,53876,54315,54741,55128,55580,56024,56448,56896,57335,57791,58278,58713,59157,59587,59978];
        ann.resp.breaths.x2 = [546,821,1301,1781,2150,2636,3058,3491,3985,4423,4861,5294,5906,6265,6771,7199,7682,8099,8532,9070,9555,10036,10532,10996,11400,11912,12408,12867,13363,13828,14356,14772,15240,15678,16106,16544,16971,17356,17752,18153,18596,19096,19539,19998,20278,20536,21038,21243,21513,22056,22518,22925,23363,23716,24139,24513,24898,25125,25526,25906,26305,26712,27155,27546,28005,28490,28923,29356,29841,30172,30694,31211,31649,32077,32467,32931,33412,33943,34465,35014,35552,35990,36407,36929,37320,37698,38125,38495,38759,39117,39534,39967,40384,40796,41228,41691,42134,42582,42989,43422,43875,44266,44688,45113,45509,45931,46417,46834,47266,47799,48237,48718,49148,49513,49919,50320,50811,51243,51755,52199,52645,53094,53479,53869,54297,54714,55141,55569,55996,56464,56891,57329,57788,58263,58712,59155,59577,59978];
    case 50
        ann.resp.breaths.x1 = [441,920,1433,1924,2407,3007,3502,4017,4504,4970,5452,5922,6396,6865,7317,7833,8341,8863,9324,9776,10228,10724,11172,11643,12126,12600,13078,13513,13996,14461,14943,15528,16037,16633,17176,17693,18189,18707,19313,19839,20287,20791,21287,21804,22365,22902,23446,24063,24593,25141,25680,26233,26743,27209,27670,28135,28648,29204,29657,30072,30576,31046,31546,31989,32446,32928,33376,33843,34413,34957,35443,36057,36600,37109,37646,38211,38828,39328,39854,40380,40946,41409,41883,42361,42839,43265,43730,44191,44687,45120,45672,46163,46680,47207,47767,48254,48680,49139,49643,50113,50587,51091,51652,52187,52680,53180,53689,54224,54772,55320,55872,56330,56857,57409,58039,58552,59083,59643];
        ann.resp.breaths.x2 = [441,910,1412,1913,2420,2989,3501,4027,4491,4982,5452,5922,6402,6850,7325,7856,8342,8880,9328,9777,10241,10743,11170,11638,12118,12609,13078,13511,13976,14466,14957,15536,16037,16639,17187,17689,18206,18734,19317,19850,20304,20784,21259,21782,22357,22909,23453,24049,24624,25141,25648,26212,26717,27208,27683,28131,28664,29181,29630,30082,30578,31063,31533,32003,32456,32916,33359,33874,34439,34951,35441,36027,36581,37093,37650,38226,38832,39313,39877,40373,40938,41416,41886,42345,42825,43274,43717,44234,44677,45129,45678,46185,46686,47214,47778,48259,48707,49164,49655,50125,50600,51122,51681,52193,52656,53178,53690,54265,54777,55315,55848,56353,56865,57466,58005,58538,59102,59656];
    case 51
        ann.resp.breaths.x1 = [80,567,985,1389,1789,2250,2689,3150,3507,3935,4317,4704,5113,5539,5983,6357,6770,7204,7602,8007,8459,8867,9233,9641,10002,10411,10854,11270,11726,12074,12452,12804,13226,13639,14100,14457,14874,15250,15646,16076,16533,16937,17324,17728,18189,18650,19096,19452,19848,20261,20704,21104,21526,21978,22348,22759,23172,23633,24015,24415,24824,25263,25689,26076,26465,26878,27257,27661,28078,28543,28996,29404,29787,30180,30524,30950,31363,31707,32080,32411,32833,33263,33702,34117,34548,34965,35409,35852,36257,36674,37052,37417,37867,38241,38693,39080,39450,39902,40267,40702,41115,41509,41900,42348,42761,43139,43583,44022,44404,44843,45233,45646,46050,46433,46837,47263,47654,48085,48493,48826,49204,49609,50000,50409];
        ann.resp.breaths.x2 = [77,583,984,1385,1781,2251,2710,3158,3507,3932,4328,4708,5119,5536,5980,6349,6771,7188,7619,8009,8463,8869,9223,9645,10030,10431,10864,11268,11728,12086,12445,12820,13226,13648,14092,14445,14857,15256,15631,16090,16533,16950,17330,17731,18190,18654,19101,19465,19850,20262,20695,21101,21528,21961,22346,22751,23168,23627,24022,24439,24830,25268,25685,26102,26464,26875,27229,27683,28100,28559,29007,29419,29794,30198,30525,30963,31369,31707,32066,32425,32831,33269,33702,34106,34555,34982,35404,35837,36265,36660,37035,37431,37867,38236,38701,39075,39466,39909,40268,40690,41091,41506,41907,42366,42772,43163,43580,44018,44419,44846,45256,45646,46047,46443,46839,47261,47662,48079,48491];
    case 52
        ann.resp.breaths.x1 = [37,450,876,1263,1689,2154,2554,2946,3354,3796,4213,4622,5052,5330,5713,6135,6565,7043,7457,7859,8254,8680,9111,9554,9898,10180,10550,10946,11370,11843,12278,12665,13070,13465,13839,14283,14713,15111,15550,15950,16333,16767,17185,17598,17985,18367,18809,19252,19626,20000,20374,20726,21113,21522,21943,22387,22802,23220,23637,24020,24420,24846,25263,25667,26059,26465,26883,27283,27683,28065,28474,28865,29296,29730,30120,30533,30941,31367,31776,32198,32624,33028,33402,33835,34243,34630,35061,35465,35891,36274,36674,37061,37478,37863,38237,38667,39093,39502,39946,40341,40772,41185,41583,41930,42261,42643,43070,43465,43939,44365,44822,45220,45602,46037,46485,46898,47320,47737,48180,48580,48970,49387,49748,50048];
        ann.resp.breaths.x2 = [40,462,900,1301,1734,2145,2541,2953,3364,3790,4222,4629,5051,5346,5726,6133,6571,7046,7452,7851,8257,8685,9112,9555,9920,10194,10558,10959,11374,11854,12276,12672,13063,13448,13844,14287,14714,15124,15551,15953,16343,16776,17187,17583,17989,18369,18821,19254,19629,20003,20357,20721,21111,21518,21945,22389,22803,23220,23637,24022,24413,24840,25268,25664,26070,26474,26875,27282,27688,28078,28469,28854,29287,29725,30119,30525,30947,31364,31781,32203,32615,33026,33406,33832,34249,34639,35061,35468,35895,36270,36666,37061,37484,37872,38231,38664,39091,39497,39946,40336,40759,41191,41585,41928,42261,42651,43063,43474,43928,44371,44820,45214,45599,46032,46491,46897,47330,47731,48174];
    case 53
        ann.resp.breaths.x1 = [220,563,893,1150,1385,1572,2007,2320,2754,3011,3254,3376,3441,3550,4035,4422,4791,5139,5509,5865,6226,6630,6965,7339,7693,8137,8493,8902,9307,9702,10098,10485,10815,11207,11639,12030,12404,12778,13122,13470,13874,14243,14643,15024,15354,15711,16102,16502,16911,17246,17615,17963,18341,18698,19039,19426,19765,20157,20500,20904,21291,21717,22057,22461,22793,23159,23511,23941,24267,24680,25080,25480,25893,26270,26622,26987,27374,27739,28143,28539,28961,29330,29743,30154,30563,30933,31337,31763,32185,32554,32876,33307,34109,34448,34865,35257,35665,36048,36478,36878,37283,37689,38054,38476,38941,39350,39746,40180,40563,41033,41513,42039,42470,42817,43257,43600,44070,44452,44865,45237,45689,46111,46489,46915];
        ann.resp.breaths.x2 = [219,530,947,1179,1380,1575,2013,2298,2773,3053,3248,3549,4053,4423,4797,5177,5510,5863,6243,6586,7009,7336,7698,8120,8505,8864,9307,9672,10078,10511,10796,11202,11659,12002,12398,12778,13121,13458,13859,14239,14635,15013,15372,15715,16111,16522,16897,17219,17615,17963,18354,18691,19064,19449,19761,20167,20473,20900,21291,21713,22061,22468,22793,23141,23521,23927,24276,24672,25094,25495,25880,26258,26627,26997,27350,27746,28142,28548,28965,29324,29735,30161,30541,30926,31348,31770,32182,32562,32910,33322,33758,34101,34460,34877,35267,35668,36053,36481,36887,37283,37687,38051,38474,38927,39355,39751,40194,40563,41038,41495,42034,42477,42820,43253,43590,44060,44424,44836,45224,45683,46111,46475,46897,47208,47646,48100,48480,48879,49270,49576,49972,50341,50742,51127,51586,52024,52415,52766,53152,53553,53954,54339,54756,55162,55584,55991,56400,56780,57197,57604,58031,58422,58849,59276,59577];
end

end

function create_readme_and_license_files(up)

fprintf('\n -- Saving ReadMe and License files')

%% License file
licence_file_path = [up.paths.data_root, 'LICENSE'];
fid = fopen(licence_file_path, 'wt');
fprintf(fid, ['The BIDMC dataset contains modified versions of PhysioNet files. The modified files are distributed under the terms for distribution of modified data files listed at: https://www.physionet.org/copying.shtml#data .',...
    '\n\n The following modifications have been made:', ...
    '\n\n - Waveform, numerics and fixed data (such as age) were downloaded from the PhysioNet MIMIC II matched waveform database. 8-min segments of data were extracted for each of the 53 ICU stays included in the BIDMC dataset.',...
    '\n\n - The data were imported into Matlab. Two annotators manually annotated breaths using the impedance signal. These manual breath annotations were included in the dataset.',...
    '\n\n - The dataset was then exported from Matlab in three formats: Matlab (r) format, in a manner which is a compatible with the RRest Toolbox of respiratory rate algorithms; (ii) CSV format; (iii) WFDB format.', ...
    '\n\n Further details of the dataset are provided in the accompanying ReadMe file.']);
fclose all;

readme_file_path = [up.paths.data_root, 'README'];
fid = fopen(readme_file_path, 'wt');

%% ReadMe file

readme_text = [...
'BIDMC Dataset',...
'\n=============',...
'\n\nTable of Contents',...
'\n=================',...
'\n\n  * [The Dataset](#the-dataset)',...
'\n  * [Revisions](#revisions)',...
'\n  * [Description](#description)',...
'\n\n# The Dataset',...
'\n\nThe BIDMC Dataset was first report in the following publication:',...
'\n\nPimentel, M.A.F. et al. Towards a Robust Estimation of Respiratory Rate from Pulse Oximeters, IEEE Transactions on Biomedical Engineering, 64(8), pp.1914-1923, 2016. [DOI: 10.1109/TBME.2016.2613124](http://doi.org/10.1109/TBME.2016.2613124) ',...
'\n\nIn this publication it was used to evaluate the performance of different algorithms for estimating respiratory rate from the pulse oximetry, or photoplethysmogram (PPG) signal. The dataset was extracted from the much larger [MIMIC II matched waveform Database](http://physionet.org/physiobank/database/mimic2wdb/matched/), which was acquired from critically-ill patients during hospital care at the Beth Israel Deaconess Medical Centre (Boston, MA, USA). The 53 recordings within the dataset, each of 8-minute duration, each contain:',...
'\n  * Physiological signals, such as the PPG, impedance respiratory signal, and electrocardiogram (ECG). These are sampled at 125 Hz.',...
'\n  * Physiological parameters, such as the heart rate (HR), respiratory rate (RR), and blood oxygen saturation level (SpO2). These are sampled at 1 Hz.',...
'\n  * Fixed parameters, such as age and gender',...
'\n  * Manual annotations of breaths. Two annotators manually annotated individual breaths in each recording using the impedance respiratory signal.',...
'\n\nThe dataset is distributed in three formats:',...
'\n  1) Matlab (r) format, in a manner which is a compatible with the [RRest Toolbox of respiratory rate algorithms](http://peterhcharlton.github.io/RRest);',...
'\n  2) CSV (comma-separated-value) format;',...
'\n  3) WFDB (WaveForm DataBase) format, which is the standard format used by [PhysioNet](https://www.physionet.org/).',...
'\n\nFor more information about the dataset, please contact the authors at:',...
'\nmarco.pimentel@eng.ox.ac.uk, peter.charlton@kcl.ac.uk .',...
'\n\n# Revisions',...
'\n\nR1: 	2017-Sept-24 		initial release',...
'\nR2:     2018-Apr-30         second release',...
'\n\n## Description',...
'\n### Matlab (r) Format',...
'\nThe *bidmc_data.mat* file contains the following subset of the dataset in a single Matlab (r) variable named *data*. The following are provided for each of the 53 recordings:',...
'\n* *ekg*:   Lead II ECG signal. Each signal is provided in a structure, where the *v* field denotes the signal values, and *fs* is the sampling frequency.',...
'\n* *ppg*:   Photoplethysmogram signal',...
'\n* *ref.resp_sig.imp*:  Impedance respiratory signal',...
'\n* *ref.breaths*:  Manual annotations of breaths provided by two independent annotators. A vector of sample numbers is provided, which correspond to the signal sample numbers.',...
'\n* *ref.params*:  Physiological parameters: *rr* (respiratory rate, derived by the monitor from the impedance signal, breaths per minute), *hr* (heart rate, derived from the ECG, beats per minute), *pr* (pulse rate, derived from the PPG, beats per minute), *spo2* (blood oxygen saturation level, %%).',...
'\n* *fix*: A structure of fixed variables, including: *id* (the MIMIC II matched waveform database subject ID and recording identifier), *loc* (the ward location), and *source* (the URLs from which the original data were downloaded).',...
'\n\n### CSV Format',...
'\n\nSeparate CSV files are provided for each recording (where ## is the subject number), containing:',...
'\n* *bidmc_##_Breaths.csv*: Manual breath annotations',...
'\n* *bidmc_##_Signals.csv*: Physiological signals',...
'\n* *bidmc_##_Numerics.csv*: Physiological parameters',...
'\n* *bidmc_##_Fix.txt*: Fixed variables',...
'\n\n### WFDB Format',...
'\n\nFive files are provided for each recording (where ## is the subject number):',...
'\n* *bidmc##.breath*: Manual breath annotations',...
'\n* *bidmc##.dat*: Waveform data file',...
'\n* *bidmc##.hea*: Waveform header file',...
'\n* *bidmc##n.dat*: Numerics data file',...
'\n* *bidmc##n.hea*: Numerics header file',...
'\n\nFurther details on the contents of each file are provided [here](https://physionet.org/tutorials/creating-records.shtml).'];

fprintf(fid, readme_text);
fclose all;

end