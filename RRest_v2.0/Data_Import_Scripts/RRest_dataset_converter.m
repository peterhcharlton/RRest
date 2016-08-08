function RRest_dataset_converter(dataset_name)
% RRest_dataset_converter converts a RRest data file from Matlab format into
% CSV and WFDB format.
%
%               RRest_dataset_converter
%
%	Inputs:
%       dataset_name    the name of the dataset for conversion, e.g.
%                       'RRsynth'.
%
%	Outputs:
%       files containing all the data written to the path
%       specified by "up.analysispath".
%
%   Requirements:
%       This requires the WFDB Toolbox, which can be downloaded from:
%           https://physionet.org/physiotools/matlab/wfdb-app-matlab/
%           
%   Further Information:
%       This version of RRest_dataset_converter is provided to convert the
%       "RRsynth" dataset originally described in:
%           Charlton P.H. and Bonnici T.B. et al. An assessment of algorithms
%           to estimate respiratory rate from the electrocardiogram and
%           photoplethysmogram, Physiological Measurement, 37(4), 2016.
%       Further information on this study can be obtained at:
%           http://peterhcharlton.github.io/RRest/yhvs_assessment.html
%       In addition, further information on RRest, including future
%       versions, can be obtained at:
%           http://peterhcharlton.github.io/RRest
%
%   Comments, Questions, Criticisms, Feedback, Contributions:
%       See: http://peterhcharlton.github.io/RRest/contributions.html
%
%   Version:
%       v.2 - published on 1st April 2016 by Peter Charlton
%               (note that this file was first uploaded on 8th August 2016)
%
%   Licence:
%       please see the accompanying file named "LICENSE"
%

fprintf('\n\n~~~~~  Starting Dataset Conversion  ~~~~~');

%% setup parameters for waveform generation
up = setup_params;

%% Load data
data = load_data(dataset_name, up);

%% Convert data to WFDB format
if up.output.wfdb_format
    convert_to_wfdb(data, dataset_name, up);
end

%% Convert data to CSV format
if up.output.csv_format
    convert_to_csv(data, dataset_name, up);
end

end

function up = setup_params

fprintf('\n - Setting up Universal Parameters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PARAMETERS TO BE SPECIFIED %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Directory from which to load the synthetic dataset:
up.analysispath = 'C:\Documents\Data\';

% direction of slash used on this computer:
up.slash = '\';

% Type(s) of output
up.output.wfdb_format = false;
up.output.csv_format = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% OTHER PARAMETERS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

up.csv.new_line = '\n';

% Check that the WFDB Matlab Toolbox has been installed
if ~exist('getWfdbClass', 'file')
    warning('Couldn''t find the WFDB Matlab Toolbox. Please install as described at the top of the file.')
end

end

function data = load_data(dataset_name, up)

fprintf('\n - Loading dataset in Matlab format');

matlab_file_path = [up.analysispath, dataset_name, 'data.mat'];

load(matlab_file_path)

end

function convert_to_wfdb(data, dataset_name, up)

fprintf('\n - Converting to WFDB format');

% lower case of dataset name to keep to usual format of PhysioNet files:
dataset_name = lower(dataset_name);
dataset_name = strrep(dataset_name, 'rrsynth', 'rrest-syn');

% Create folder to save data in
data_type = 'wfdb';
save_folder = [up.analysispath, dataset_name, '_', data_type, up.slash];
if ~exist(save_folder, 'dir')
    mkdir(save_folder)
end
cd(save_folder)

% add licence file
licence_file_path = [save_folder, 'LICENSE'];
fid = fopen(licence_file_path, 'wt');
fprintf(fid, 'This dataset is licensed under a Creative Commons Attribution 4.0 International License. Further details are available at: https://creativecommons.org/licenses/by/4.0/');
fclose all;

% cycle through each subject
no_subjs = length(data);
file_names = cell(no_subjs,1);
for subj_no = 1 : no_subjs
    
    % setup command
    data_mat(:,1) = data(subj_no).ppg.v; data_mat(:,2) = data(subj_no).ekg.v*1000;
    file_names{subj_no} = [dataset_name, sprintf('%.3d', subj_no)];
    fs = data(subj_no).ekg.fs;
    units = 'au/mV';
    descrip = ['<age>:  <sex>:  <diagnoses>: (none)  <medications>: (none)  <ventilation>: ' data(subj_no).fix.ventilation, '  <subject group>: ' data(subj_no).group  '  <heart rate>: ' num2str(data(subj_no).ref.params.hr.v), ' ', data(subj_no).ref.params.hr.units '  <respiratory rate>: ' num2str(data(subj_no).ref.params.rr.v), ' ', data(subj_no).ref.params.rr.units];
    sig_names = {['Pleth, ' data(subj_no).ppg.method], ['II, ' data(subj_no).ekg.method]};
    
    % convert this subject's data into WFDB format
    mat2wfdb(data_mat, file_names{subj_no}, fs, [], units, descrip, [], sig_names);
    
    clear data_mat fs units descrip sig_names
    
end

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

% create a DBS file
file_name = 'DBS';
file_path = [save_folder, file_name];
fid = fopen(file_path, 'w');
fprintf(fid, [dataset_name, '\t', 'Synthetic ECG and PPG, modulated by respiration']);

end

function convert_to_csv(data, dataset_name, up)

fprintf('\n - Converting to CSV format');

% lower case of dataset name to keep to usual format of PhysioNet files:
dataset_name = lower(dataset_name);
dataset_name = strrep(dataset_name, 'rrsynth', 'rrest-syn');

% Create folder to save data in
data_type = 'csv';
save_folder = [up.analysispath, dataset_name, '_', data_type, up.slash];
if ~exist(save_folder, 'dir')
    mkdir(save_folder)
end
cd(save_folder)

% add licence file
licence_file_path = [save_folder, 'LICENSE'];
fid = fopen(licence_file_path, 'wt');
fprintf(fid, 'This dataset is licensed under a Creative Commons Attribution 4.0 International License. Further details are available at: https://creativecommons.org/licenses/by/4.0/');
fclose all;

% cycle through each subject
no_subjs = length(data);
file_names = cell(no_subjs,1);
for subj_no = 1 : no_subjs
    
    % setup command
    data_mat(:,1) = data(subj_no).ppg.v; data_mat(:,2) = data(subj_no).ekg.v*1000;
    file_names{subj_no} = [dataset_name, sprintf('%.3d', subj_no)];
    fs = data(subj_no).ekg.fs;
    units = 'au/mV';
    descrip = ['<age>:  <sex>:  <diagnoses>: (none)  <medications>: (none)  <ventilation>: ' data(subj_no).fix.ventilation, '  <subject group>: ' data(subj_no).group  '  <heart rate>: ' num2str(data(subj_no).ref.params.hr.v), ' ', data(subj_no).ref.params.hr.units '  <respiratory rate>: ' num2str(data(subj_no).ref.params.rr.v), ' ', data(subj_no).ref.params.rr.units];
    sig_names = {['Pleth, ' data(subj_no).ppg.method], ['II, ' data(subj_no).ekg.method]};
    
    % create a file of fixed variables for this subject
    curr_filename = [file_names{subj_no}, '_fix.txt'];
    fid = fopen(curr_filename, 'wt');
    text_to_write = [file_names{subj_no}, up.csv.new_line, ...
        'Signals: ' sig_names{1} '; and ' sig_names{2}, up.csv.new_line, ...
        'Sampling frequency: ' num2str(fs) ' Hz', up.csv.new_line, ...
        'Respiratory modulation type: ' upper(data(subj_no).group), up.csv.new_line, ...
        'Simulated respiratory rate: ' num2str(data(subj_no).ref.params.rr.v), ' ', data(subj_no).ref.params.rr.units, up.csv.new_line, ...
        'Simulated heart rate: ' num2str(data(subj_no).ref.params.hr.v), ' ', data(subj_no).ref.params.hr.units];
    fprintf(fid, text_to_write); clear text_to_write
    
    % create a file of data for this subject
    curr_filename = [file_names{subj_no}, '_data.csv'];
    csvwrite(curr_filename, data_mat); clear data_mat
    
    % close file
    fclose(fid); clear fid
    
    clear fs units descrip sig_names
    
end

end