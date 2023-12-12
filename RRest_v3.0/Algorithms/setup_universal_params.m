function up = setup_universal_params(period_orig)
% SETUP_UNIVERSAL_PARAMS sets the universal parameters for use throughout
% "RRest".
%
%               setup_universal_params('vortal_factors')
%
%	Inputs:
%		period_orig - the name of the dataset on which to do the
%		analysis, e.g. if the file is called 'vortal_factors_data.mat',
%		then the dataset's name is 'vortal_factors'
%
%	Outputs:
%       up      - a struct of universal parameters
%
%   Context:    called by "RRest.m"
%           
%   Further Information:
%       This version is specifically designed to facilitate reproduction of
%       the analysis performed in:
%           Charlton P.H. et al. Extraction of respiratory signals from the 
%           electrocardiogram and photoplethysmogram: technical and physiological
%           determinants, Physiological Measurement, 38(5), 2017
%           DOI: https://doi.org/10.1088/1361-6579/aa670e
%       Further information on this study can be obtained at:
%           http://peterhcharlton.github.io/RRest/factors_assessment.html
%       In addition, further information on RRest, including future
%       versions, can be obtained at:
%           http://peterhcharlton.github.io/RRest/index.html
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% FILE PATHS TO BE SPECIFIED %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% (requires editing) %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Please note that your system may require the slashes in file paths to be
% of the opposite direction. In which case, change the following:
up.paths.slash_direction = filesep;     % usually a backslash for Windows, forward slash for Linux

% Specify path of data root folder. For instance, if you specify
% "up.paths.root_folder = 'C:\Documents\Data\';", then data will be saved
% in the directory located at 'C:\Documents\Data\DATASETNAME\', where
% "DATASETNAME" is the name of the dataset being analysed, e.g. "vortal_rest".
up.paths.root_folder = 'C:\DataBase\BIDMC\bidmc-ppg-and-respiration-dataset-1.0.0\';%'C:\Documents\Data\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% LOAD PREVIOUS UNIVERSAL PARAMETERS %%%%%%%%
%%%%%%%%%%%%%%%% (no editing required) %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Provide licence details
provide_licence_details
% set current directory to that of this file
cd(fileparts(mfilename('fullpath')))
% Add all functions corresponding to options for each component to the path
addpath(genpath(fileparts(mfilename('fullpath'))));
% Specify name of universal params file:
up_filename = 'up';
up.paths.root_data_folder = [up.paths.root_folder, period_orig, up.paths.slash_direction];
up.paths.data_save_folder = [up.paths.root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
up_path = [up.paths.data_save_folder, up_filename];
% Specify whether or not to re-make the universal_params (needs to be
% remade if any settings have changed)
redo_up = true;
% ------------------------------
% Don't remake the universal_params file unless needed:
if exist([up_path, '.mat'], 'file') && ~redo_up
    fprintf('\n--- Loading Universal Parameters ');
    load(up_path);
    return
end
% Otherwise make the universal_params file:
fprintf('\n--- Creating Universal Parameters ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% ALGORITHMS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% (edit if you wish) %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the stages of the algorithms (best left alone):
up.al.key_components = {'extract_resp_sig', 'estimate_rr', 'fuse_rr'};      % To run the analysis in full this should be: {'extract_resp_sig', 'estimate_rr', 'fuse_rr'}
% Specify methods for extraction of respiratory signals (feature / filter, ecg / ppg):
up.al.options.extract_resp_sig = {'ppg_feat','ekg_feat'};%ppg_filt', 'ekg_feat', 'ekg_filt'};  % Possible methods: 'ppg_feat', 'ekg_feat', 'ppg_filt', 'ekg_filt'
% ------------------------------------------------------
% Specify the components for feature-based extraction of respiratory signals:
up.al.sub_components.ppg_feat = {'EHF', 'PDt', 'FPt', 'FMe', 'RS', 'ELF'};  % Should read: {'EHF', 'PDt', 'FPt', 'FMe', 'RS', 'ELF'}
up.al.sub_components.ekg_feat = {'EHF', 'RDt', 'FPt', 'FMe', 'RS', 'ELF'};  % Should read: {'EHF', 'RDt', 'FPt', 'FMe', 'RS', 'ELF'}
% Specify the interchangeable technique(s) to be used for each component of feature-based extraction of respiratory signals:
up.al.options.PDt = {'IMS'};                                                % Possible methods: 'DCl', 'COr', 'IMS'
up.al.options.RDt = {'GC'};                                                 % Possible methods: 'GC', 'ME'
up.al.options.FMe = {'am', 'fm', 'bw', 'pyppg_bm','Xb200'};                                     % Possible methods: 'pulW', 'am', 'fm', 'bw', 'pk', 'on', 'bwm', 'qrsW', 'qrsA', 'pca', 'qrS', 'rsS', 'Rang'
% ------------------------------------------------------
up.al.options.RS = {'linB'};                                                % Possible methods: 'cub', 'cubB', 'brg', 'lin', 'brgB', 'linB'
% Specify the interchangeable technique(s) to be used for filter-based respiratory signal extraction:
up.al.options.ekg_filt = {'Wfm', 'Wam', 'BFi'};                             % Possible methods: 'Wfm', 'Wam', 'CCF', 'BFi'
up.al.options.ppg_filt = {'Wfm', 'Wam', 'BFi'};                             % Possible methods: 'Wfm', 'Wam', 'CCF', 'BFi'
% Specify the interchangeable technique(s) for RR Estimation
%up.al.options.estimate_rr = {'FTS', 'ARS', 'ARM', 'ACF', 'WCH', 'PKS', 'ZeX', 'PZX', 'CtO', 'CtA'};       % Possible methods: 'FTS', 'ARS', 'ARM', 'ARP', 'ARPz', 'ACF', 'WCH', 'PKS', 'ZeX', 'PZX', 'CtO', 'CtA'
up.al.options.estimate_rr = {'CtO'};   
% Different methods for fusion of RR estimates:
up.al.options.fuse_rr = {'fus_mod', 'fus_temp'};                            % Possible methods: 'fus_mod', 'fus_temp'
% Components for each method of extraction of RR fusion:
%up.al.sub_components.fus_mod = {'SFu', 'SPA'};                              % Possible methods: 'SFu', 'PMC', 'PRC', 'SPA'
up.al.sub_components.fus_mod = {'SFu'};  
%up.al.sub_components.fus_temp = {'TFu'};                                    % Possible methods: 'TFu'
up.al.sub_components.fus_temp = {};   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% DIRECTORIES %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% (no editing required) %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify file paths
up.paths.data_load_filename = [period_orig, '_data'];
up.paths.data_load_filename2 = [period_orig, 'data'];
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
up.paths.filenames.cc = '_cc';
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
up.paths.filenames.comb_results_table = 'results_table';
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

abs_path=pwd;
end_p=strfind(abs_path,'Peter_toolbox')-1;
abs_path=abs_path(1:end_p);
bm_dir=[abs_path,'python\PPG_feats\all_BM\'];
up.paths.bm = bm_dir;

fpt_dir=[abs_path,'python\PPG_feats\Fiducial_points\'];
up.paths.fpt = fpt_dir;

% extract names of biomarkers
if sum(strcmp(up.al.options.FMe, 'pyppg_bm'))
    dataset_name = up.paths.data_load_filename;
    dataset_name = strrep(dataset_name, 'data', '');
    num_id=sprintf('%02d', 1);
    loadpath = [up.paths.bm, dataset_name, num_id,'.mat'];
    load(loadpath);
    bm_names = fieldnames(all_bm);
    bm_names = bm_names(:)';
%     bm_names=[];%bm_names(24);
    bm_names(strcmp(bm_names, 'Index of pulse')) = [];
    bm_names(strcmp(bm_names, 'TimeStamp')) = [];
    up.al.options.FMe = [up.al.options.FMe, bm_names];
    up.al.options.FMe(strcmp(up.al.options.FMe, 'pyppg_bm')) = [];
    up.al.options.bm_names = bm_names;
end

% Make directories
directories_to_make = {up.paths.data_load_folder, up.paths.data_save_folder, up.paths.results_folder, up.paths.plots_save_folder, up.paths.tables_save_folder};
for s = 1 : length(directories_to_make)
    rel_path = directories_to_make{s};
    if ~exist(rel_path, 'dir')
        mkdir(rel_path)
    end    
end
% Move data for processing into the relevant folder (you may wish to delete the original file)
data_orig_path = [up.paths.root_folder, up.paths.data_load_filename, '.mat'];
data_new_path = [up.paths.data_load_folder, up.paths.data_load_filename, '.mat'];
data_orig_path2 = [up.paths.root_folder, up.paths.data_load_filename2, '.mat'];
if exist(data_orig_path, 'file')
    copyfile(data_orig_path, data_new_path);
elseif exist(data_orig_path2, 'file')
    copyfile(data_orig_path2, data_new_path);
end

% path of the folder in which to store results figures. This can be
% wherever you like.
up.paths.paper_figures_folder = up.paths.plots_save_folder;   % in which to save eps figures for direct import into publication.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% ANALYSIS PARAMETERS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% (no editing required) %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find out what subjects need processing
data_path = [up.paths.data_load_folder, up.paths.data_load_filename, '.mat'];
up.analysis.run_analysis = true;
if exist(data_path, 'file')
    load(data_path);
else
    warning('Note there was no raw data for this period...')
    path00 = [up.paths.root_folder, 'vortal_rest', up.paths.slash_direction];
    path0 = [path00, 'Analysis_files', up.paths.slash_direction, 'Data_for_Analysis', up.paths.slash_direction];
    path2 = 'vortal_rest_data';
    load([path0, path2]);
    up.analysis.run_analysis = false;
    up.paths.equipment_type = '';                                           % Possible equipment types: either '_clin' for clinical monitor, or empty, '', for raw signal acquisition
end
up.paramSet.subj_list = 1:length(data);%:5

% Find out which sub and sub-sub group analyses need conducting
if sum(strcmp(fieldnames(data),'group'))
    up.paramSet.groups = extractfield(data, 'group');
    if sum(strcmp(fieldnames(data), 'sub_group'))
        up.paramSet.sub_groups = extractfield(data, 'sub_group');
    end
end

% Find out what ppg and ekg signals there are:
possible_sigs = fieldnames(data);
up.paramSet.ppg_sigs = possible_sigs(~cellfun(@isempty, strfind(possible_sigs, 'ppg')));
up.paramSet.ekg_sigs = possible_sigs(~cellfun(@isempty, strfind(possible_sigs, 'ekg')));
up.paramSet.scg_sigs = possible_sigs(~cellfun(@isempty, strfind(possible_sigs, 'scg')));

% window parameters
up.paramSet.winLeng = 32;                                                   % the duration of each window in secs
up.paramSet.winStep = 0;                                                    % the number of secs between each consecutive window
up.paramSet.buffer_period = 10;                                             % the number of secs to ignore at the start of each resp sig (since the filters might not have stabilised).

% Filter characteristics: Eliminate VHFs (above frequency content of signals)
up.paramSet.elim_vhf.ppg.Fpass = 38.5;  % in HZ
up.paramSet.elim_vhf.ppg.Fstop = 33.12;  % in HZ   (33.12 and 38.5 provide a -3 dB cutoff of 35 Hz)
up.paramSet.elim_vhf.ppg.Dpass = 0.05;
up.paramSet.elim_vhf.ppg.Dstop = 0.01;
up.paramSet.elim_vhf.ekg.Fpass = 103.7;  % in HZ
up.paramSet.elim_vhf.ekg.Fstop = 98;  % in HZ    (98 and 103.7 provide a -3 dB cutoff of 100 Hz)
up.paramSet.elim_vhf.ekg.Dpass = 0.05;
up.paramSet.elim_vhf.ekg.Dstop = 0.01;
% Filter characteristics: Eliminate Mains frequencies (50 Hz)
up.paramSet.elim_mains.fcuts = [41, 47.2, 52.8, 59];  % in HZ (41, 47.2, 52.8, 59 provide a -3 dB cutoff of 45-55 Hz)
up.paramSet.elim_mains.Dpass = 0.05;
up.paramSet.elim_mains.Dstop = 0.01;
% Filter characteristics: Eliminate HFs (above resp freqs)
up.paramSet.elim_hf.Fpass = 1.2;  % in Hz
up.paramSet.elim_hf.Fstop = 0.8899;  % in Hz     (1.2 and 0.8899 provide a -3dB cutoff of 1 Hz)
up.paramSet.elim_hf.Dpass = 0.05;
up.paramSet.elim_hf.Dstop = 0.01;
% up.paramSet.elim_hf.Fpass = 0.661;  % in Hz
% up.paramSet.elim_hf.Fstop = 0.567;  % in Hz     (0.661 and 0.567 provide a -3dB cutoff of 0.6 Hz)
% up.paramSet.elim_hf.Dpass = 0.05;
% up.paramSet.elim_hf.Dstop = 0.01;
% Filter characteristics: Eliminate HFs (above resp freqs) - for ref RR estimation
up.paramSet.elim_hf_ref_rr.Fpass = 1.2;  % in Hz
up.paramSet.elim_hf_ref_rr.Fstop = 0.8899;  % in Hz     (1.2 and 0.8899 provide a -3dB cutoff of 1 Hz)
up.paramSet.elim_hf_ref_rr.Dpass = 0.05;
up.paramSet.elim_hf_ref_rr.Dstop = 0.01;
% Filter characteristics: Eliminate LFs (below cardiac freqs): For 30bpm cutoff
up.paramSet.elim_sub_cardiac.Fpass = 0.63;  % in Hz
up.paramSet.elim_sub_cardiac.Fstop = 0.43;  % in Hz     (0.63 and 0.43 provide a - 3dB cutoff of 0.5 Hz)
up.paramSet.elim_sub_cardiac.Dpass = 0.05;
up.paramSet.elim_sub_cardiac.Dstop = 0.01;
% Filter characteristics: Eliminate VLFs (below resp freqs): For 4bpm cutoff
up.paramSet.elim_vlf.Fpass = 0.157;  % in Hz
up.paramSet.elim_vlf.Fstop = 0.02;   % in Hz     (0.157 and 0.02 provide a - 3dB cutoff of 0.0665 Hz)
up.paramSet.elim_vlf.Dpass = 0.05;
up.paramSet.elim_vlf.Dstop = 0.01;
% Filter characteristics: duration of Tukey window taper in secs
up.paramSet.tukey_win_duration_taper = 2;

% frequency at which to resample feature-based respiratory signals
up.paramSet.resample_fs = 5;
% frequency at which to sample filter-based respiratory signals
up.paramSet.filt_resample_fs = 25;

% Plausible range of HRs [bpm]
up.paramSet.hr_range = [30, 220];

% Plausible range of plausible RRs [bpm]
up.paramSet.rr_range = [4 60];

% AR model order
up.paramSet.ar_model_order = 8;
up.paramSet.ar_resample_freq = 5;

% FFT downsample freq
up.paramSet.fft_resample_freq = 5;

% decide whether or not to redo stats
up.analysis.redo_stats = 1;

% QRS detection windows in secs
up.paramSet.qrs_detect.win_length = 10;
up.paramSet.qrs_detect.win_overlap = 5;

% PPG Peak Detector specifications
% CO
smooth_filt = load('parafilt'); % Load PPG smoothing filter 
up.paramSet.CO_peak_det.PPG_smoother.coeffs = smooth_filt.parafiltPPG;
up.paramSet.CO_peak_det.PPG_smoother.details = smooth_filt.details;
up.paramSet.CO_peak_det.upctl = 0.9;
up.paramSet.CO_peak_det.lpctl = 0.1;

% Determine how to obtain reference RRs:
ref_fields = fieldnames(data(1).ref);
if sum(strcmp(ref_fields, 'resp_sig')) && sum(strcmp(fieldnames(data(1).ref.resp_sig), 'paw'))
    up.paramSet.ref_method = 'paw';
elseif sum(strcmp(ref_fields, 'resp_sig')) && sum(strcmp(fieldnames(data(1).ref.resp_sig), 'imp'))
    up.paramSet.ref_method = 'imp_sqi';  % using the impedance SQI
    % up.paramSet.ref_method = 'imp_agree'; % using the agreement between RRs derived from impedance using time- and freq-domain methods
elseif sum(strcmp(ref_fields, 'resp_sig')) && sum(strcmp(fieldnames(data(1).ref.resp_sig), 'unknown'))
    up.paramSet.ref_method = 'imp_sqi';  % using the impedance SQI
elseif sum(strcmp(ref_fields, 'resp_sig')) && sum(strcmp(fieldnames(data(1).ref.resp_sig), 'band'))
    up.paramSet.ref_method = 'band_agree'; % using the agreement between RRs derived from impedance using time- and freq-domain methods
elseif sum(strcmp(ref_fields, 'breaths'))
    up.paramSet.ref_method = 'breaths';
elseif sum(strcmp(ref_fields, 'params')) && sum(strcmp(fieldnames(data(1).ref.params), 'rr'))
    up.paramSet.ref_method = 'rrs';
end

% Determine whether or not to calculate CCs:
if sum(strcmp(up.paramSet.ref_method, 'rrs')) || sum(strcmp(up.paramSet.ref_method, 'breaths'))
    up.paramSet.calc_CCs = false;
else
    up.paramSet.calc_CCs = true;
end

% Determine whether to calculate stats for impedance:
if sum(strcmp(fieldnames(data(1).ref), 'params')) && sum(strcmp(fieldnames(data(1).ref.params), 'rr')) ...
        && sum(strcmp(fieldnames(data(1).ref.params.rr), 'method')) && ~isempty(strfind(data(1).ref.params.rr.method, 'impedance'))
    up.analysis.imp_stats = true;
else
    up.analysis.imp_stats = false;
end
% Determine whether to calculate stats for synthetic data:
if ~isempty(strfind(up.paths.root_data_folder, 'synth'))
    up.analysis.calc_synth_stats = true;
else
    up.analysis.calc_synth_stats = false;
end
clear data

% Ref RR threshold
up.paramSet.resp_sig_thresh.paw = 0.42;                                     % thresh for use of paw (just a single paw, not 2*paw)
up.paramSet.resp_filter_roll_off = 1+(1/3);

% Paw SNR
up.paramSet.paw_snr_thresh = -8;

% Reference RR method
up.paramSet.ref_rr_method = 'paw_thresh';                                   % For Vortal: sum_both_thresh, imp_thresh, paw_thresh; For MIMIC: comb_imp; For CAPNOBASE: co2_zex

% Impedance RR methods
up.paramSet.imp_rr_method = 'num';   % num or thresh

% Threshold for Q- and S- points (either side of R)
up.paramSet.q_s_thresh = 0.1;  % in secs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% SAVE UNIVERSAL PARAMETERS %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% (no editing required) %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(up_path, 'up');

end

function provide_licence_details

licence_details = ['\n\n RRest  Copyright (C) 2017  King''s College London',...
    '\n This program comes with ABSOLUTELY NO WARRANTY', ... 
    '\n and is available under the GNU public license.', ...
    '\n For details see the accompanying LICENSE.txt file.\n\n'];

fprintf(licence_details)

end