function run_vortal_downsampler
% RUN_VORTAL_DOWNSAMPLER provides extra signals for the 'vortal_young_elderly'
% dataset with a range of sampling frequencies.
%
% This script is required to produce the results pertaining to sampling
% frequency reported in the following publication:
%
%       Charlton P.H. et al. Extraction of respiratory signals from the 
%       electrocardiogram and photoplethysmogram: technical and physiological
%       determinants, Physiological Measurement, 38(5), 2017
%       DOI: https://doi.org/10.1088/1361-6579/aa670e
%
%	Inputs:
%		vortal_young_elderly_data.mat - the relevant dataset, which is available
%		at:  http://peterhcharlton.github.io/RRest/vortal_dataset.html
%       DOI: %%%%%%%%%%  TBC  %%%%%%%%%%%%%%
%
%	Outputs:
%       vortal_factors_data.mat - the 'vortal_factors' dataset, which is
%       simply a copy of the 'vortal_young_elderly' dataset, with the addition 
%       of some downsampled signals.
%           
%   Further Information:
%       This function is part of v.3 of the RRest toolbox. It is provided
%       to facilitate reproduction of the analysis performed in the above
%       publication. Further information on this study can be obtained at:
%           http://peterhcharlton.github.io/RRest/factors_assessment.html
%
%       In addition, further information on RRest, including future
%       versions, can be obtained at:
%           http://peterhcharlton.github.io/RRest
%
%   Comments, Questions, Criticisms, Feedback, Contributions:
%       See: http://peterhcharlton.github.io/RRest/contributions.html
%
%   Version:
%       v.3 - published on 4th May 2017 by Peter H. Charlton
%
%   Licence:
%       Available under the GNU public license - please see the accompanying
%       file named "LICENSE"
%

%% Provide licence details
provide_licence_details

%% Setup Universal Parameters
up = setup_universal_params;

%% Load vortal_ext_data file
data = load_data_file(up);

%% Downsample PPG signal
data = downsample_signal(data, 'ppg', up);

%% Downsample ECG signal
data = downsample_signal(data, 'ekg', up);

%% Save new data file
save_data(data, up);

end

function provide_licence_details

% Print licence details at the command window:

licence_details = ['\n\n run_vortal_downsampler', ...
    '\n Copyright (C) 2017  King''s College London',...
    '\n This program comes with ABSOLUTELY NO WARRANTY', ... 
    '\n and is available under the GNU public license.', ...
    '\n For details see the accompanying LICENSE.txt file.\n\n'];

fprintf(licence_details)

end

function up = setup_universal_params

fprintf('\n--- Setting up universal parameters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% IDENTIFY DATA FILE %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% (no editing required) %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[data_file_name,up.paths.root_data_folder,~] = uigetfile('','Please select the ''vortal_young_elderly_data.mat'' file');
up.paths.data_filepath = [up.paths.root_data_folder, data_file_name];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% SET UNIVERSAL PARAMETERS %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% (no editing required) %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

new_dataset_name = 'vortal_factors_data';
up.paths.data_save_filepath = [up.paths.root_data_folder, new_dataset_name];

% Specify sampling frequencies to be downsampled to
up.pub.freqs.ppg = 125./[2,4,6,8,10,15];
up.pub.freqs.ekg = 500./[2,3,4,5,7,10];                                     % may also require 25 Hz to replicate analyses in paper

% Specify signals to be downsampled (those from the clinical monitor)
up.pub.rel_sigs.ppg = 'ppgclin';
up.pub.rel_sigs.ekg = 'ekgclin';

end

function data = load_data_file(up)

fprintf('\n--- Loading data');

% load data file
load(up.paths.data_filepath)

end

function data = downsample_signal(data, sig_type, up)

fprintf(['\n--- Downsampling ' sig_type ' signals']);

% for each subject
for subj_no = 1 : length(data)
    
    % identify relevant signal
    eval(['rel_sig = up.pub.rel_sigs.' sig_type ';'])
    eval(['orig_sig = data(subj_no).' rel_sig ';'])
    orig_sig.t = [0:(length(orig_sig.v)-1)]./orig_sig.fs;
    
    % cycle through the required sampling frequencies
    eval(['rel_freqs = up.pub.freqs.' sig_type ';'])
    for fs = rel_freqs
        
        % downsample
        downsample_factor = round(orig_sig.fs/fs);
        new_sig.v = downsample(orig_sig.v, downsample_factor);
        new_sig.t = downsample(orig_sig.t, downsample_factor);
        new_sig.fs = fs;
        
        % re-interpolate back up to original sampling frequency
        new_sig.v = interp1(new_sig.t, new_sig.v, orig_sig.t(orig_sig.t<=new_sig.t(end)), 'spline');
        new_sig.t = orig_sig.t(orig_sig.t<=new_sig.t(end));
        new_sig.fs = orig_sig.fs;
        
        % store downsampled signal
        sig_name = [rel_sig, 'DS', num2str(round(fs)) 'Hz'];
        eval(['data(subj_no).' sig_name ' = new_sig;'])
        
    end
end

end

function save_data(data, up)

fprintf('\n--- Saving data');

% save data
save(up.paths.data_save_filepath, 'data', '-v7.3')

end
