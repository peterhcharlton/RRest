function RRSYNTH_data_generater
% RRSYNTH_data_generater generates synthetic ECG and PPG data under the
% influence of each of three respiratory modulations: baseline wander (BW),
% amplitude modulation (AM), and frequency modulation (FM).
%
%               RRSYNTH_data_generater
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
%   Further Information:
%       This version of RRSYNTH_data_generater is provided to facilitate
%       reproduction of the analysis performed in:
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
%               (note that this file was first uploaded on 30th June 2016)
%
%   Licence:
%       please see the accompanying file named "LICENSE"
%

fprintf('\n\n~~~~~  Starting Synthetic Data Generation  ~~~~~');

%% setup parameters for waveform generation
up = setup_params;

%% Load data
[ecg_data, ppg_data] = load_sample_beat_data;

%% Generate range of data samples covering different mods and respiratory rates
[ecg_data, ppg_data] = generate_synthetic_data(ecg_data, ppg_data, up);

%% Save modulated data in common format
export_to_common_format(ecg_data, ppg_data, up);

end

function up = setup_params

fprintf('\n - Setting up Universal Parameters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PARAMETERS TO BE SPECIFIED %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Directory in which to save the simulated dataset:
up.analysispath = 'C:\Documents\Data\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   OTHER PARAMETERS   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   (don't need editing)   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Duration of simulated data:
up.no_secs = 210;

% Amplitude of data:
up.sig_range = 1;

% Sampling frequency:
up.fs = 500;

% Magnitude of modulations:
up.a_r = 0.1;           % (for BW and AM)
up.a_r_fm = up.a_r/2;   % (for FM)

% Specify range of RRs and HRs for which signals should be generated:
up.range_plausible_rrs = [4, 60];   % e.g. RRs from 4 - 60 bpm
up.rr_interval = 2;                 %       at an interval of 2 bpm
up.range_plausible_hrs = [30, 200]; % e.g. HRs from 30 - 200 bpm
up.hr_interval = 5;                 %       at an interval of 5 bpm

% Specify the constant HRs and RRs at which to simulate the data
up.const_hr = 80;
up.const_rr = 20;

end

function [ecg_data, ppg_data] = load_sample_beat_data

fprintf('\n - Loading sample beat data');

% ECG beat data
ecg_data.beat.t = 0:0.01:1;
ecg_data.beat.v = [-0.120000000000000,-0.120000000000000,-0.125000000000000,-0.122500113257034,-0.125000000000000,-0.120000000000000,-0.125000000000000,-0.117500065569862,-0.125000000000000,-0.117499958278698,-0.120000000000000,-0.120000000000000,-0.124999856955537,-0.122500005960186,-0.105000131124091,-0.115000000000000,-0.110000000000000,-0.122499946358326,-0.125000000000000,-0.109999558946239,-0.100000357611158,-0.0899996065808297,-0.0750000000000000,-0.0775001728453928,-0.0850000000000000,-0.0975001490224131,-0.100000000000000,-0.0925001251639052,-0.0949997377518178,-0.0824998986647591,-0.0800000000000000,-0.0800000000000000,-0.100000429133389,-0.0299993562231754,0.464997210632973,0.997500745023246,0.954992990821313,0.134999821173104,-0.194999845035165,-0.165000035761116,-0.119999880796281,-0.117499958273724,-0.110000202646322,-0.107499934437955,-0.109999928477768,-0.100000000000000,-0.100000000000000,-0.0900000000000000,-0.0949999761592562,-0.0824998390749790,-0.0800000000000000,-0.0675001609250210,-0.0649999761592562,-0.0574998867564668,-0.0500000000000000,-0.0400000000000000,-0.0300000000000000,-0.0124999344379545,0.00500020264632248,0.0275000417262755,0.0649997615925615,0.100000035761115,0.145000309929670,0.189999988078207,0.224999666229586,0.252499970199070,0.270000000000000,0.275000000000000,0.250000214566695,0.197500232447252,0.114999761592563,0.0325003040057223,-0.0299992132554528,-0.0849997496721898,-0.105000023840744,-0.130000000000000,-0.130000000000000,-0.140000000000000,-0.130000000000000,-0.127500196709585,-0.120000000000000,-0.120000000000000,-0.110000000000000,-0.110000000000000,-0.105000000000000,-0.105000000000000,-0.100000000000000,-0.100000000000000,-0.100000000000000,-0.100000000000000,-0.104999821194421,-0.110000000000000,-0.110000000000000,-0.105000000000000,-0.114999773512934,-0.120000000000000,-0.115000000000000,-0.122500113257034,-0.115000000000000,-0.120000000000000,-0.110000000000000];

% PPG beat data
ppg_data.beat.t = 0:0.01:1;
ppg_data.beat.v = [1077,1095.27676998271,1150.39696200268,1253.64011384634,1410.01430649103,1614.23738042021,1852.88862847947,2107.19963850395,2354.10827203910,2571.83504395768,2747.58603743220,2878.78047565119,2967.12445239160,3018.75510162723,3042.94721680823,3047.45001341082,3037.62026345592,3018.14128598677,2991.80738548651,2959.68011205818,2924.12507822977,2885.14270668177,2843.86782834153,2802.08903767062,2759.59973773619,2716.75011175682,2673.89981522322,2630.12661808970,2584.21488466353,2533.75439619717,2480.71246051141,2423.10076799285,2360.92007271860,2294.87110683952,2226.11520653275,2155.41275815825,2085.84032425344,2019.82560572195,1959.08788519998,1905.61966740180,1860.15018328118,1823.59858526554,1799.16523439130,1785.68124422590,1781.40496155451,1785.65000894081,1795.60758813888,1810.23123770638,1825.88000834451,1843.33369762174,1863.68739383102,1877.70370775466,1882.41001847768,1880.67124988824,1874.06744888836,1862.89381761287,1847.14991953269,1829.89768909817,1810.34993145378,1788.72101296379,1765.94996721702,1742.03613261809,1717.39751028193,1692.75877484727,1668.68005245276,1644.48098438338,1619.84239547012,1595.97728854980,1572.41003606020,1548.84233116767,1525.27517210550,1501.70741230815,1478.27021279132,1455.77367020325,1433.54980777827,1412.12496274662,1390.69993145582,1369.27500149014,1348.40755237670,1328.05378822197,1308.39978542052,1289.88870123678,1272.51984025750,1256.29379973178,1241.20988734577,1229.73332273703,1217.97235166145,1206.26252309710,1195.54992698555,1185.27003397509,1176.98743219884,1169.12989658769,1160.98994575908,1153.99245954403,1146.99245931931,1138.20624497094,1127.13495738213,1112.99225070036,1097.13740724184,1083.85480956071,1080];

end

function [ecg_data, ppg_data] = generate_synthetic_data(ecg_data, ppg_data, up)

fprintf('\n - Generating synthetic data');

% Designed to generate ECG and PPG signals for the following:
% - constant HR = 80, with RR varying from 4 to 60 at 2 bpm intervals;
% - constant RR = 20, with HR varying from 30 to 200 at 5 bpm intervals;
% - each modulation on its own (BW, AM and FM)

% Setup
mods = {'bw', 'am', 'fm'};
sigs = {'ppg_data', 'ecg_data'};

% Cycle through each signal
for sig_no = 1 : length(sigs)
    seg_counter = 0;
        
    % Holding each of HR and RR constant in turn
    for const = {'hr', 'rr'}
                    
        %% settings
        % find HR and RR frequencies to be used
        if strcmp(const, 'hr')
            settings.hr.f_bpm = up.const_hr;
            settings.rr.f_bpm = up.range_plausible_rrs(1) : up.rr_interval : up.range_plausible_rrs(2);
        else
            settings.hr.f_bpm = up.range_plausible_hrs(1) : up.hr_interval : up.range_plausible_hrs(2);
            settings.rr.f_bpm = up.const_rr;
        end
        settings.rr.f_hz = settings.rr.f_bpm./60;
        settings.hr.f_hz = settings.hr.f_bpm./60;
        settings.rr.w = 2*pi*settings.rr.f_hz;
        settings.hr.w = 2*pi*settings.hr.f_hz;
        
        % cycle through mods
        for mod_no = 1 : length(mods)
            
            % cycle through HRs
            for hr_no = 1 : length(settings.hr.w)
                
                %% Find un-modulated data for this HR
                temp.hr = settings.hr.f_bpm(hr_no);
                [ppg_data, ecg_data] = find_unmod_data(ppg_data, ecg_data, temp.hr, up);
                
                % Cycle through RRs
                for rr_no = 1 : length(settings.rr.w)
                    
                    %% Find modulated data for this combination of settings
                    % Find out what the combination of settings is:
                    temp.sig = sigs{sig_no};
                    temp.mod = mods{mod_no};
                    temp.rr = settings.rr.w(rr_no);
                    seg_counter = seg_counter+1;
                    
                    eval(['unmod_data = ' temp.sig '.unmod;']);
                    mod_data = generate_modulated_data(unmod_data, temp, up);
                    
                    sig_data{seg_counter,1} = mod_data.v;
                    sig_time{seg_counter,1} = mod_data.t;
                    sig_f_r(seg_counter,1) = settings.rr.f_hz(rr_no);
                    sig_f_c(seg_counter,1) = settings.hr.f_hz(hr_no);
                    sig_mod{seg_counter,1} = mods{mod_no};
                    
                end
            end
        end
        
        % Store the simulated signals
        eval([temp.sig, '.mod.v = sig_data;']);
        eval([temp.sig, '.mod.t = sig_time;']);
        eval([temp.sig, '.mod.f_r = sig_f_r;']);
        eval([temp.sig, '.mod.f_c = sig_f_c;']);
        eval([temp.sig, '.mod.type = sig_mod;']);
    end
end
end

function [ppg_data, ecg_data] = find_unmod_data(ppg_data, ecg_data, hr, up)

%% Find un-modulated data for this HR
for sig = {'ppg_data', 'ecg_data'}
    % Extract data
    eval(['beat_data = ' sig{1,1} '.beat;']);
    % Normalise data
    beat_data.v = up.sig_range * (beat_data.v - min(beat_data.v))/(range(beat_data.v));
    % Find duration of a single beat
    beat_duration = 60/hr;
    beat_samples = beat_duration*up.fs;
    % resample data at appropriate sampling rate
    int = (1/(beat_samples-1));
    t_old = 0:(1/100):1;
    t_new = 0:int:1;
    beat_data.v = interp1(t_old(1:(end-1)), beat_data.v(1:(end-1)), t_new, 'pchip');
    beat_data.t = t_new;
    % Repeat data
    tot_no_samples = up.no_secs*up.fs;
    no_reps = ceil(tot_no_samples/length(beat_data.t));
    beat_data.seg.v = repmat(beat_data.v, [1,no_reps]);
    beat_data.seg.v = beat_data.seg.v(1:tot_no_samples);
    beat_samples = length(beat_data.t);    % incase the number of samples isn't a whole number
    beat_data.seg.t = [0:(tot_no_samples-1)]*(1/beat_samples)*(beat_samples/up.fs);
    % Re-insert data
    eval([sig{1,1} '.unmod = beat_data.seg;']);
end

end

function mod_data = generate_modulated_data(unmod_data, temp, up)

% setup
w_r = temp.rr;
mod = temp.mod;

% modulate the unmodulated signal using the specified modulation
switch mod
    case 'bw'
        mod_data.v = unmod_data.v + up.a_r*sin(w_r*unmod_data.t);
        mod_data.t = unmod_data.t;
    case 'am'
        mod_data.v = detrend(unmod_data.v).*((1/up.a_r)+cos(w_r*unmod_data.t));
        mod_data.t = unmod_data.t;
    case 'fm'
        temp.t = unmod_data.t + up.a_r_fm*sin(w_r*unmod_data.t);
        temp.v = unmod_data.v;
        mod_data.t = unmod_data.t;
        mod_data.v = interp1(temp.t, temp.v, mod_data.t);
end

%% Eliminate nans
mod_data.t = mod_data.t(~isnan(mod_data.t));
mod_data.v = mod_data.v(~isnan(mod_data.v));

end

function export_to_common_format(ecg_data, ppg_data, up)

fprintf('\n - Saving data');

%% Create a data structure based on the Karlen data format

for subj_el = 1:length(ecg_data.mod.v)
    
    % setup
    if exist('data', 'var')
        struct_el = length(data)+1;
    else
        struct_el = 1;
    end
    
    % Insert fixed params
    data(1,struct_el).group = ppg_data.mod.type{subj_el};
    data(1,struct_el).fix.ventilation = 'simulated';
    data(1,struct_el).ref.params.hr.v = 60*ppg_data.mod.f_c(subj_el);
    data(1,struct_el).ref.params.hr.units = 'beats/min';
    data(1,struct_el).ref.params.hr.method = 'simulated constant HR throughout recording';
    data(1,struct_el).ref.params.rr.v = 60*ppg_data.mod.f_r(subj_el);
    data(1,struct_el).ref.params.rr.units = 'breaths/min';
    data(1,struct_el).ref.params.rr.method = 'simulated constant RR throughout recording';
    
    % insert PPG signal
    data(1,struct_el).ppg.v = ppg_data.mod.v{subj_el};
    data(1,struct_el).ppg.fs = up.fs;
    data(1,struct_el).ppg.method = 'simulated finger PPG';
    
    % insert EKG signal
    data(1,struct_el).ekg.v = ecg_data.mod.v{subj_el};
    data(1,struct_el).ekg.fs = up.fs;
    data(1,struct_el).ekg.method = 'simulated ECG';
    
    % Generate breath timings
    mod_data_t = [0:(length(data(subj_el).ppg.v)-1)]*(1/(up.fs));
    breaths.t = mod_data_t(1):(1/ppg_data.mod.f_r(subj_el)):mod_data_t(end);
    
    % insert reference RRs
    data(1,struct_el).ref.breaths.t = breaths.t;
    data(1,struct_el).ref.breaths.method = 'simulated';
    data(1,struct_el).ref.breaths.units.t = 's';
    
end

% Save to file
save([up.analysispath, 'RRSYNTHdata'], 'data')

end