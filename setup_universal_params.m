function up = setup_universal_params(period_orig)
%SETUP_UNIVERSAL_PARAMS sets up the parameters for use throughout the function
%"run_rr_algorithms.m".
%
%               setup_universal_params
%
%	Inputs:
%		period_orig - the period of Vortal recording on which to do the analysis, e.g. 'rest', 'walk', 'ex', 'rec', 'synth'
%
%	Outputs:
%       up      - a struct of universal parameters
%
%   See also:
%       run_rr_algorithms.m
%

%% Computer specific path settings
% Change these if you are using this for the first time (they're at the
% bottom of this file)
up = computer_specific_settings;

%% Compulsory setup tasks

% set current directory to that of this file
cd(fileparts(mfilename('fullpath')))
% Add all functions corresponding to options for each component to the path
addpath(genpath(fileparts(mfilename('fullpath'))));
% Specify name of universal params file:
up_foldername = 'universal_params';
up_filename = 'up';
up_path = [fileparts(mfilename('fullpath')), up.paths.slash_direction, up_foldername, up.paths.slash_direction, up_filename];
% Specify whether or not to re-make the universal_params
redo_up =1;

%% Optional setup tasks

% Don't remake the universal_params file unless needed:
if exist([up_filename, '.mat'], 'file') & ~redo_up
    fprintf('\n--- Loading Universal Parameters ');
    load(up_filename);
    return
end

% Otherwise make the universal_params file:
fprintf('\n--- Creating Universal Parameters ');

%% Specify the components of RR algorithms to be tested
% This specifies the structure of algorithms (best left alone):
up.al.key_components = {'extract_resp_sig', 'estimate_rr', 'fuse_rr'}; % {'extract_resp_sig', 'estimate_rr', 'fuse_rr'};  	%   % To run the analysis in full this should be: {'extract_resp_sig', 'estimate_rr', 'fuse_rr'}
% Different methods for extraction of respiratory signals (feature / filter, ecg / ppg):
up.al.options.extract_resp_sig = {'ekg_feat', 'ppg_feat', 'ekg_filt', 'ppg_filt'}; % To run the analysis in full this should be: {'ppg_feat', 'ekg_feat', 'ekg_filt', 'ppg_filt'}
up.al.options.extract_resp_sig = {'ppg_feat', 'ekg_feat'};
% Components for each method of extraction of respiratory signals:
up.al.sub_components.ppg_feat = {'EHF', 'PDt', 'FPt', 'FMe', 'RS', 'ELF'};  % Should read: {'EHF', 'PDt', 'FPt', 'FMe', 'RS', 'ELF'}
up.al.sub_components.ekg_feat = {'EHF', 'RDt', 'FPt', 'FMe', 'RS', 'ELF'};  % Should read: {'EHF', 'RDt', 'FPt', 'FMe', 'RS', 'ELF'}
% Interchangeable techniques for each component:
up.al.options.PDt = {'IMS'};                                                % Possible methods: 'DCl', 'COr', 'IMS'
up.al.options.RDt = {'GC'};                                                 % Possible methods: 'GC', 'ME'
up.al.options.FMe = {'pulW', 'am', 'fm', 'bw', 'pk', 'on', 'bwm', 'qrsW', 'qrsA', 'pca'};  % Possible methods: 'pulW', 'am', 'fm', 'bw', 'pk', 'on', 'bwm', 'qrsW', 'qrsA', 'pca'
up.al.options.FMe = {'am', 'fm', 'bw'};
up.al.options.RS = {'linB'};                                                % Possible methods: 'cub', 'cubB', 'brg', 'lin', 'brgB', 'linB'
up.al.options.ekg_filt = {'Wfm', 'Wam', 'CCF', 'BFi'};               % Possible methods: 'Wfm', 'Wam', 'CCF', 'BFi'
up.al.options.ppg_filt = {'Wfm', 'Wam', 'CCF', 'BFi'};               % Possible methods: 'Wfm', 'Wam', 'CCF', 'BFi'
up.al.options.estimate_rr = {'FTS', 'ACF', 'ARPz', 'WCH', 'ARM', 'ARP', 'CtO', 'CtA', 'ARS', 'PKS', 'ZeX', 'PZX'}; %  Possible methods: 'ARM', 'ACF', 'ARP', 'ARPz', 'CtO', 'CtA', 'ARS', 'FTS', 'PKS', 'ZeX', 'PZX'
up.al.options.estimate_rr = {'FTS', 'CtO'}; % , 'ARe'
% Different methods for fusion of RR estimates:
up.al.options.fuse_rr = {'fus_mod', 'fus_temp'};                            % Possible methods: 'fus_mod', 'fus_temp'
%up.al.options.fuse_rr = {'fus_mod'};
% Components for each method of extraction of RR fusion:
up.al.sub_components.fus_mod = {'SPA', 'SFu', 'PMC', 'PRC'};                       % Possible methods: 'SFu', 'PMC', 'PRC'
%up.al.sub_components.fus_mod = {'SFu'};
up.al.sub_components.fus_temp = {'TFu'};                                    % Possible methods: 'TFu'

%% Specify the folders in which to save results files

% Synthetic Data
period = period_orig;
if strcmp(period, 'synth')
    up.paths.root_data_folder = [up.paths.root_folder, 'RRSYNTH', up.paths.slash_direction];
    up.paths.data_load_filename = 'RRSYNTHdata';
elseif strcmp(period, 'mimic')
    % MIMIC Data
    up.paths.root_data_folder = [up.paths.root_folder, 'MIMIC_MATCHED', up.paths.slash_direction];
    up.paths.data_load_filename = 'MIMICIIdata';
elseif strcmp(period, 'capnobase')
    % CapnoBase Data
    up.paths.root_data_folder = [up.paths.root_folder, 'CAPNOBASE', up.paths.slash_direction];
    up.paths.data_load_filename = 'CAPNOBASEdata';
elseif strcmp(period, 'listen')
    % LISTEN data
    up.paths.root_data_folder = [up.paths.root_folder, 'RRLISTEN', up.paths.slash_direction];
    up.paths.data_load_filename = 'RRLISTENdata';
else
    % VORTAL Data
    up.paths.equipment_type = '';                                           % Possible equipment types: either '_clin' for clinical monitor, or empty, '', for raw signal acquisition
    if strcmp(period, 'rest')
        up.paths.root_data_folder = [up.paths.root_folder, 'VORTAL', up.paths.equipment_type, up.paths.slash_direction];
    else
        eval(['up.paths.root_data_folder = [up.paths.root_folder, ''VORTAL', up.paths.equipment_type, '_' upper(period), up.paths.slash_direction, '''];']);
    end
    eval(['up.paths.data_load_filename = ''VORTAL_' lower(period) '_data', up.paths.equipment_type, ''';']);    
end

% File paths
up.paths.data_load_folder = [up.paths.root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Data_for_Analysis', up.paths.slash_direction];
up.paths.data_save_folder = [up.paths.root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
up.paths.results_folder = [up.paths.root_data_folder, 'Analysis_files' , up.paths.slash_direction, 'Results', up.paths.slash_direction];
up.paths.plots_save_folder = [up.paths.results_folder, 'Figures', up.paths.slash_direction];
up.paths.tables_save_folder = [up.paths.results_folder, 'Tables', up.paths.slash_direction];
up.paths.filenames.respSigs = '_respSigs.mat';
up.paths.filenames.int_respSigs = '_int_respSigs.mat';
up.paths.filenames.pulse_peaks = '_PDt';
up.paths.filenames.qrss = '_RDt';
up.paths.filenames.fid_pts = '_FPt';
up.paths.filenames.elim_vhf = '_EHF';
up.paths.filenames.feat_meas = '_FMe';
up.paths.filenames.resampler = '_RS';
up.paths.filenames.elim_vlf2 = '_ELF';
up.paths.filenames.filt = '_flt';
up.paths.filenames.rrEsts = '_rrEsts';
up.paths.filenames.fuse = 'fuse';
up.paths.filenames.ref_rrs = '_rrRef';
up.paths.filenames.sqi = '_sqi';
up.paths.filenames.stats = '_stats';
up.paths.filenames.win_timings = '_wins';
up.paths.filenames.global_stats = 'stats_global';
up.paths.filenames.ehv_stats = 'stats_ehv';
up.paths.filenames.yhv_stats = 'stats_yhv';
up.paths.filenames.group_stats = 'stats_';
up.paths.filenames.global_data = 'data_global';
up.paths.filenames.yhv_data = 'data_yhv';
up.paths.filenames.ehv_data = 'data_ehv';
up.paths.filenames.group_data = 'data_';
up.paths.filenames.global_BA = 'BA_global';
up.paths.filenames.imp_BA = 'BA_imp';
up.paths.filenames.yhv_BA = 'BA_yhv';
up.paths.filenames.ehv_BA = 'BA_ehv';
up.paths.filenames.group_BA = 'BA_';
up.paths.filenames.study_stats = 'study_stats';
up.paths.filenames.alg_names = 'alg_names';
up.paths.filenames.imp_alg_names = 'imp_alg_names';
up.paths.filenames.results_table = 'results_table';
up.paths.filenames.BA_results_table = 'BA_results_table';
up.paths.filenames.stats_results_table = 'stats_results_table';
up.paths.filenames.subj_results_table = 'subj_results_table';
up.paths.filenames.win_results_table = 'win_results_table';
up.paths.filenames.win_results_headers = 'win_results_headers';
up.paths.filenames.win_data = 'win_data';
up.paths.filenames.win_data_imp = 'win_data_imp';
up.paths.filenames.win_data_raw_clin = 'win_data_raw_clin';
up.paths.filenames.win_data_raw_clin_imp = 'win_data_raw_clin_imp';
up.paths.filenames.vortal_results = 'vortal_results';
up.paths.filenames.synth_results = 'synth_results';
up.paths.filenames.hr_rr_scatter = 'hr_rr_scatter';
up.paths.filenames.stacked_bar = 'stacked_bar';
up.paths.filenames.joint_vortal_BA = 'BA_rest_and_rec';
up.paths.filenames.rest_and_rec_win_data = 'rest_and_rec_win_data';
up.paths.filenames.temp_prec_plot ='temp_prec_plot';
up.paths.filenames.feat_filt_plot = 'feat_filt_plot';

% Make directories
directories_to_make = {up.paths.data_load_folder, up.paths.data_save_folder, up.paths.results_folder, up.paths.plots_save_folder, up.paths.tables_save_folder};
for s = 1 : length(directories_to_make)
    rel_path = directories_to_make{s};
    if ~exist(rel_path, 'dir')
        mkdir(rel_path)
    end    
end
% Move data for processing into the relevant folder
data_orig_path = [up.paths.root_data_folder, up.paths.data_load_filename, '.mat'];
data_new_path = [up.paths.data_load_folder, up.paths.data_load_filename, '.mat'];
if exist(data_orig_path, 'file')
    copyfile(data_orig_path, data_new_path);
end
   % you may wish to delete the original file

%% Specify analysis parameters

% Find out what subjects need processing
data_path = [up.paths.data_load_folder, up.paths.data_load_filename, '.mat'];
if exist(data_path, 'file')
    load(data_path);
else
    warning('Note there was no raw data for this period...')
    path00 = [up.paths.root_folder, 'VORTAL', up.paths.slash_direction];
    path0 = [path00, 'Analysis_files', up.paths.slash_direction, 'Data_for_Analysis', up.paths.slash_direction];
    eval(['path2 = ''VORTAL_rest_data'';']);
    load([path0, path2]);
end
up.paramSet.subj_list = 1:length(data);
up.paramSet.groups = extractfield(data, 'group');
clear data

% window parameters
up.paramSet.winLeng = 32;                                                   % the duration of each window in secs
up.paramSet.winStep = 0;                                                    % the number of secs between each consecutive window
%up.paramSet.refRR_winStep = 3;                                              % the number of secs between each consecutive window over which the ref RR is calculated
up.paramSet.buffer_period = 10;                                             % the number of secs to ignore at the start of each resp sig (since the filters might not have stabilised).

% Filter characteristics
% Eliminate VHFs (above frequency content of signals)
up.paramSet.elim_vhf.ppg.Fpass = 38.5;  % in HZ
up.paramSet.elim_vhf.ppg.Fstop = 33.12;  % in HZ   (33.12 and 38.5 provide a -3 dB cutoff of 35 Hz)
up.paramSet.elim_vhf.ppg.Dpass = 0.05;
up.paramSet.elim_vhf.ppg.Dstop = 0.01;
up.paramSet.elim_vhf.ekg.Fpass = 103.7;  % in HZ
up.paramSet.elim_vhf.ekg.Fstop = 98;  % in HZ    (98 and 103.7 provide a -3 dB cutoff of 100 Hz)
up.paramSet.elim_vhf.ekg.Dpass = 0.05;
up.paramSet.elim_vhf.ekg.Dstop = 0.01;
% Eliminate Mains frequencies (50 Hz)
up.paramSet.elim_mains.fcuts = [41, 47.2, 52.8, 59];  % in HZ (41, 47.2, 52.8, 59 provide a -3 dB cutoff of 45-55 Hz)
up.paramSet.elim_mains.Dpass = 0.05;
up.paramSet.elim_mains.Dstop = 0.01;

% Eliminate HFs (above resp freqs)
up.paramSet.elim_hf.Fpass = 1.2;  % in Hz
up.paramSet.elim_hf.Fstop = 0.8899;  % in Hz     (1.2 and 0.8899 provide a -3dB cutoff of 1 Hz)
up.paramSet.elim_hf.Dpass = 0.05;
up.paramSet.elim_hf.Dstop = 0.01;
% Eliminate LFs (below cardiac freqs): For 30bpm cutoff
up.paramSet.elim_sub_cardiac.Fpass = 0.63;  % in Hz
up.paramSet.elim_sub_cardiac.Fstop = 0.43;  % in Hz     (0.63 and 0.43 provide a - 3dB cutoff of 0.5 Hz)
up.paramSet.elim_sub_cardiac.Dpass = 0.05;
up.paramSet.elim_sub_cardiac.Dstop = 0.01;
% Eliminate VLFs (below resp freqs): For 4bpm cutoff
up.paramSet.elim_vlf.Fpass = 0.157;  % in Hz
up.paramSet.elim_vlf.Fstop = 0.02;   % in Hz     (0.157 and 0.02 provide a - 3dB cutoff of 0.0665 Hz)
up.paramSet.elim_vlf.Dpass = 0.05;
up.paramSet.elim_vlf.Dstop = 0.01;
% duration of Tukey window taper in secs
up.paramSet.tukey_win_duration_taper = 2;
% QRS detection wins in secs
up.paramSet.qrs_detect.win_length = 10;
up.paramSet.qrs_detect.win_overlap = 5;

% frequency at which to resample feature-based respiratory signals
up.paramSet.resample_fs = 5;
% frequency at which to sample filter-based respiratory signals
up.paramSet.filt_resample_fs = 25;

% CWT parameters
up.paramSet.hr_range = [30, 220];

% RR estimation specifications
% Range of plausible RR freqs
up.paramSet.rr_range = [4 60];
% AR model order
up.paramSet.ar_model_order = 8;
up.paramSet.ar_resample_freq = 5;
% FFT downsample freq
up.paramSet.fft_resample_freq = 5;
% Ref RR threshold
up.paramSet.resp_sig_thresh.sum_both = 0.60; % thresh for use of paw and imp combined
up.paramSet.resp_sig_thresh.imp = 0.1; % thresh for use of imp  (just a single imp, not 2*imp)
up.paramSet.resp_sig_thresh.paw = 0.42; % thresh for use of paw (just a single paw, not 2*paw)
up.paramSet.resp_filter_roll_off = 1+(1/3);
% Paw SNR
up.paramSet.paw_snr_thresh = -8;
% Reference RR method
up.paramSet.ref_rr_method = 'paw_thresh'; % For Vortal: sum_both_thresh, imp_thresh, paw_thresh; For MIMIC: comb_imp; For CAPNOBASE: co2_zex
%up.paramSet.ref_rr_method = 'timings';
up.paramSet.ref_rr_method = 'comb_imp';
% Impedance RR methods
up.paramSet.imp_rr_method = 'num';   % num or thresh

% decide whether or not to redo stats
up.analysis.redo_stats = 1;

% PPG Peak Detector specifications
% CO
smooth_filt = load('parafilt'); % Load PPG smoothing filter 
up.paramSet.CO_peak_det.PPG_smoother.coeffs = smooth_filt.parafiltPPG;
up.paramSet.CO_peak_det.PPG_smoother.details = smooth_filt.details;
up.paramSet.CO_peak_det.upctl = 0.9;
up.paramSet.CO_peak_det.lpctl = 0.1;

% DClifton
opts.delta_fraction = 0.75;
opts.delta_median_window = 5;
opts.delta_decay_time = 10;     % Was originally on 100, changed to include early peaks in time series (as the first few mins of the Paw signal weren't being picked up).
opts.delta_decay_factor = 0.5;
opts.delta_minimum = 0.2;
opts.init_window_size = 50;
opts.init_buffer_size = 1800000;
opts.delta_decay_time = 50;
opts.delta_fraction = 1;
up.paramSet.DClifton_opts = opts; clear opts

%% save universal params file
save(up_path, 'up');

end

function up = computer_specific_settings

%% These are the settings which need to be changed when running this on a new computer

% Specify path of data root folder
up.paths.slash_direction = '\';
up.paths.root_folder = 'C:\Documents\Data\';

% specify other folders
up.paths.paper_figures_folder = 'C:\Users\pc13\Dropbox\VORTAL\VORTAL_theoret_lims_yhvs\Figures\';   % in which to save eps figures for direct import into publication.
up.paths.db_data = 'C:\Documents\Data\VORTAL\Analysis_files\Processed_Data\db_data.mat';

% % Add additional methods (non-published) path for Pete's computer
% cname = getenv('computername');
% if strcmp(cname, 'BIOENG043-PC')
%     addpath('C:\Documents\Matlab analysis\AR_RR\');
%     up.al.options.ekg_filt = {'ARa', 'ARf', 'Wfm', 'Wam'}; %, 'ARb'};
%     up.al.options.ppg_filt = {'ARa', 'ARf', 'Wfm', 'Wam'}; %, 'ARb', , 'BFi', 'CCF'
%     addpath(genpath('C:\Users\pc13\Documents\GitHub\phd\Attractor_Reconstruction\'));
% end

end