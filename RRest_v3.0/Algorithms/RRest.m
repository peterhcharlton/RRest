function RRest(period)
% RRest runs RR algorithms on ECG and PPG signals using each possible
% combination of options, as specified in "setup_universal_params.m".
%
%               RRest('vortal_rest')
%
%	Inputs:
%		data            data files should be stored in the specified format
%                       in the directory specified in
%                       "setup_universal_params.m". Data can be downloaded
%                       using the "
%       period          this string specifies the dataset to be analysed.
%                       Only the 'mimic' dataset has been used with this
%                       version of the toolbox.
%
%	Outputs:
%       for each subject, N, the following files are made:
%           N_int_respSigs      intermediate respiratory signals,
%           N_respSigs          final respiratory signals
%           N_rrEsts            RR estimates
%           N_rrRef             Reference RR values
%           N_sqi               Signal Quality Index values
%       for the entire dataset, the following files are made:
%           alg_names           Names of RR algorithms tested
%           win_data            Data for every algorithm tested, every
%                               window, and every subject.
%
%   Context:    This is the main file used to run the algorithms. It calls
%               lots of other functions, contained in separate files.
%           
%   Further Information:
%       This version of the RRest is provided to facilitate reproduction of
%       the analysis performed in:
%           Charlton P.H. et al. Extraction of respiratory signals from the 
%           electrocardiogram and photoplethysmogram: technical and physiological
%           determinants, Physiological Measurement, 38(5), 2017
%           DOI: https://doi.org/10.1088/1361-6579/aa670e
%       Further information on this study can be obtained at:
%           http://peterhcharlton.github.io/RRest/factors_assessment.html
%       In addition, further information on RRest, including future
%       versions, can be obtained at:
%           http://peterhcharlton.github.io/RRest
%
%   Comments, Questions, Criticisms, Feedback, Contributions:
%       See: http://peterhcharlton.github.io/RRest/contributions.html
%
%   Version:
%       v.3 - published on 4th May 2017 by Peter Charlton
%
%   Licence:
%       please see the accompanying file named "LICENSE"
%

%% Setup Universal Parameters
% The universal parameters are used throughout the algorithms
up = setup_universal_params(period);

up.analysis.run_analysis = 1;
if up.analysis.run_analysis
    
    %% Estimate RRs from ECG and PPG
    % Carry out processing for each stage of the algorithms
    for key_comp_no = 1 : length(up.al.key_components)
        feval(up.al.key_components{key_comp_no}, up);
    end
    
    %% Conduct Signal Quality Assessment of Signals
    % Each window of ECG and PPG is quality assessed
    calculate_sqi(up);
    
    %% Estimate Reference RRs
    estimate_ref_rr(up);
    
    %% Calculate resp signal correlation coefficients
    if up.paramSet.calc_CCs
        calculate_resp_sig_cc(up);
    end
    
end

%% Statistical Analysis
calc_stats(up);
create_table_of_algorithms_comb(up);

end