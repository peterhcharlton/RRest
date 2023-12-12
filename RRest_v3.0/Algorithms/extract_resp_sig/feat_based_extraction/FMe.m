function FMe(option, up)
%FMe measures features from peak and trough values as specified in
% PC's literature review.
%	            FMe(option, up)
%
%	Inputs:
%		option          the option which has led to this function being used
%       up              universal parameters structure
%
%	Outputs:
%       ...
%

fprintf('\n--- Measuring Features ');
log_int_respSig = 1;             % Has value 1 unless this is a final respiratory signal

for subj = up.paramSet.subj_list
    
    sig_type = option(1:3);
    %% Cycle through each signal of this type
    eval(['sigs = up.paramSet.' sig_type '_sigs;']);
    for sig_no = 1 : length(sigs)
        curr_sig = sigs{sig_no};
        
        %% Skip if this processing has been done previously
        iden_resp_sig_file_ending
        savepath = [up.paths.data_save_folder, num2str(subj), ending];
        filecontents = whos('-file', savepath);
        var_names = extractfield(filecontents, 'name');
        rel_log = zeros(size(var_names));
        for s = 1 : length(var_names)
            if strfind(var_names{s}, [curr_sig, up.paths.filenames.fid_pts])
                rel_log(s) = 1;
            end
        end
        rel_var_names = var_names(logical(rel_log));
        for rel_var_name_no = 1 : length(rel_var_names)
            % could move biomarker loading to here
            for current_opt_no = 1 : length(up.al.options.FMe)
                % skip out ones that aren't relevant to the ppg
                if strcmp(sig_type, 'ppg') && ( strcmp(up.al.options.FMe{current_opt_no}, 'qrsW') ...
                        || strcmp(up.al.options.FMe{current_opt_no}, 'qrsA') ...
                        || strcmp(up.al.options.FMe{current_opt_no}, 'qrS') ...
                        || strcmp(up.al.options.FMe{current_opt_no}, 'rsS') ...
                        || strcmp(up.al.options.FMe{current_opt_no}, 'Rang'))
                    continue
                end
                
                % skip out ones that aren't relevant to the ekg
                if strcmp(sig_type, 'ekg') && strcmp(up.al.options.FMe{current_opt_no}, 'pulW')
                    continue
                end
                
                temp = strfind(rel_var_names{rel_var_name_no},'_');
                start_el = temp(1); clear temp
                eval(['save_name = ''' curr_sig, up.paths.filenames.feat_meas, up.al.options.FMe{current_opt_no}, rel_var_names{rel_var_name_no}(start_el:end) ''';']); clear start_el
                exist_log = check_exists(savepath, save_name);
                if exist_log
                    continue
                end
                
                %% Load relevant data
                
                % load file containing biomarkers from pyPPG
                dataset_name = up.paths.data_load_filename;
                dataset_name = strrep(dataset_name, 'data', '');
                num_id=sprintf('%02d', subj);
                loadpath = [up.paths.bm, dataset_name, num_id,'.mat'];
                load(loadpath);

                loadpath = [up.paths.fpt, dataset_name, num_id,'*.mat'];
                file = dir(loadpath);
                loadpath = [up.paths.fpt, file.name];
                load(loadpath);
                    
                % only load signal data if not loaded already
                if ~exist('sig_data', 'var')
                    loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.int_respSigs];
                    rel_name = [curr_sig, up.paths.filenames.elim_vhf];
                    load(loadpath, rel_name);
                    eval(['sig_data = ' rel_name ';']);
                end
                
                % only load fiducial point data if not loaded already
                if ~exist('rel_data', 'var')
                    loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.int_respSigs];
                    rel_name = rel_var_names{rel_var_name_no};
                    load(loadpath, rel_name);
                    eval(['rel_data = ' rel_name ';']);
                    %% Settings
                    % choose which peaks and onsets:
                    rel_data.tr = rel_data.tr_min;
                    rel_data.p = rel_data.p_max;
                    
                    if isempty(rel_data.tr.t) || isempty(rel_data.p.t)
                        [peaks.t, peaks.v, onsets.t, onsets.v] = deal([]);
                    else
                        
                        % Want an onset followed by a peak.
                        if rel_data.tr.t(1) > rel_data.p.t(1)
                            rel_data.p.t = rel_data.p.t(2:end);
                            rel_data.p.v = rel_data.p.v(2:end);
                        end
                        % Want the same number of peaks and onsets
                        diff_in_length = length(rel_data.p.t) - length(rel_data.tr.t);
                        if diff_in_length > 0
                            rel_data.p.t = rel_data.p.t(1:(end-diff_in_length));
                            rel_data.p.v = rel_data.p.v(1:(end-diff_in_length));
                        elseif diff_in_length < 0
                            rel_data.tr.t = rel_data.tr.t(1:(end-diff_in_length));
                            rel_data.tr.v = rel_data.tr.v(1:(end-diff_in_length));
                        end
                        % find onsets and peaks
                        onsets.t = rel_data.tr.t; onsets.t = onsets.t(:);
                        onsets.v = rel_data.tr.v; onsets.v = onsets.v(:);
                        peaks.t = rel_data.p.t; peaks.t = peaks.t(:);
                        peaks.v = rel_data.p.v; peaks.v = peaks.v(:);
                        
                        % exclude ectopics
                        % following: Mateo, J. & Laguna, P., 2003. Analysis of heart rate variability in the presence of ectopic beats using the heart timing signal. IEEE transactions on bio-medical engineering, 50(3), pp.334–43. Available at: http://www.ncbi.nlm.nih.gov/pubmed/12669990.
                        if ~isempty(strfind(up.paths.data_load_filename, 'vortal'))
                            tk_neg1 = peaks.t(1:(end-2));
                            tk = peaks.t(2:(end-1));
                            tk_pos1 = peaks.t(3:end);
                            r = 2*abs( (tk_neg1 - (2*tk) + tk_pos1)./ ...
                                ( (tk_neg1-tk).*(tk_neg1 - tk_pos1).*(tk-tk_pos1) ) );
                            thresh = min([4.3*std(r), 0.5]);
                            
                            %%%%%%%%%%%%%%%%%%%%%%
                            % additional rule inserted by PC:
                            % thresh = 0.5;   % so that artificial data with a very low variability doesn't trigger.
                            %%%%%%%%%%%%%%%%%%%%%%
                            temp = [0;r;0]; temp = logical(temp>thresh);
                            
                            tk_neg1 = onsets.t(1:(end-2));
                            tk = onsets.t(2:(end-1));
                            tk_pos1 = onsets.t(3:end);
                            r = 2*abs( (tk_neg1 - (2*tk) + tk_pos1)./ ...
                                ( (tk_neg1-tk).*(tk_neg1 - tk_pos1).*(tk-tk_pos1) ) );
                            thresh = min([4.3*std(r), 0.5]);
                            %%%%%%%%%%%%%%%%%%%%%%
                            % additional rule inserted by PC:
                            %thresh = 0.5;   % so that artificial data with a very low variability doesn't trigger.
                            %%%%%%%%%%%%%%%%%%%%%%
                            temp2 = [0;r;0]; temp2 = logical(temp2>thresh);
                            
                            peaks.t(temp | temp2) = nan; peaks.v(temp | temp2) = nan;
                            onsets.t(temp | temp2) = nan; onsets.v(temp | temp2) = nan;
                            %                     % to check:
                            %                     old_peaks  = rel_data.p; old_onsets = rel_data.tr;
                            %                     plot(diff(old_peaks.t), 'b'), hold on, plot(diff(peaks.t), 'r')
                            %                     close all
                            %                     plot(diff(old_onsets.t)), hold on, plot(diff(onsets.t))
                            %                     close all
                        end
                        clear temp temp2 r tk tk_neg1 tk_pos1 thresh
                    end
                end
                
                %% Measure Features
                sig_data.wave_type = curr_sig;   % for PCA

                tmp_name=up.al.options.FMe{current_opt_no};
                if sum(strcmp(up.al.options.bm_names, tmp_name))
                    % Define current biomarker                             
                    bm_name=up.al.options.FMe{current_opt_no};
                    bm=eval(['[all_bm(:).',bm_name,'];'])';

                    % Define peaks and onsets
                    fs=sig_data.fs;
                    mode_on='pym';
                    mode_sp='pym';
                    [on, sp]=get_on_and_sp(PPG_fiducials,onsets,peaks,mode_on,mode_sp,sig_data);

                    % Find rigth fiducial indexes
                    ind_sp=find(~isnan(sp.t) & ~isnan(sp.v) & ~isnan(sp.i));%find(~isnan(sp));
                    ind_on=find(~isnan(sp.t) & ~isnan(sp.v) & ~isnan(sp.i));%find(~isnan(on));
                    fp_ind=sp.i(intersect(ind_sp,ind_on));
                    bm_ind=find(~isnan(bm));
                    [~, bmi, fpi]= intersect(bm_ind,fp_ind);
                    for i=1:3
                        fns=cell2mat(fieldnames(sp));
                        fn=fns(i);
                        new_sp.(fn)=sp.(fn)(fpi);
                        new_on.(fn)=on.(fn)(fpi);
                    end

                    % Find rigth biomarker Time Stamps
                    bm=bm(bmi);
                    TimeStamp=[all_bm.('TimeStamp')]';
                    TimeStamp=TimeStamp(bm_ind);

                    % Define final biomarkers
                    BM.t = mean([new_on.t, new_sp.t], 2);
                    BM.v = bm;
                    
                    % Normalise
                    feat_data.v = double(BM.v./nanmean(BM.v));
                    feat_data.t = double(BM.t);
                    feat_data.fs = fs;
                elseif tmp_name=="Xb200"
                    mode_on='pyc';
                    mode_sp='pym';
                    [new_on, new_sp]=get_on_and_sp(PPG_fiducials,onsets,peaks,mode_on,mode_sp,sig_data);
                    feat_data = feval(['calc_bw'], new_sp, new_on, sig_data.fs, sig_data, up);
                else
                    feat_data = feval(['calc_' tmp_name], peaks, onsets, sig_data.fs, sig_data, up);
                end
                feat_data.timings = rel_data.timings;
                
                %% Save results
                eval([save_name ' = feat_data;']);
                save_or_append_data
            end
            clear rel_data
        end
        clear sig_data
    end
end

end

%% Fuction to get new onsets and peaks
function [new_on, new_sp]=get_on_and_sp(PPG_fiducials,onsets,peaks,mode_on,mode_sp,sig_data)
    fs=sig_data.fs;
    sp_py=[];
    on_py=[];
    sp_py.t=(double([PPG_fiducials.sp]')+1)/fs;
    sp_py.v=sig_data.v([PPG_fiducials.sp]+1);
    sp_py.i=[PPG_fiducials.('Index of pulse')]';
    on_py.t=(double([PPG_fiducials.on]')+1)/fs;
    on_py.v=sig_data.v([PPG_fiducials.on]+1);
    on_py.i=[PPG_fiducials.('Index of pulse')]';

    on_pyc=on_py;
    for ti=2:length(on_py.t)-1
        tempon=find(onsets.t>on_py.t(ti-1) & onsets.t<on_py.t(ti+1));
        if length(tempon)>0
            [~,tind]=min(abs(on_py.t(ti)-onsets.t(tempon)));
            on_pyc.t(ti,1)=onsets.t(tempon(tind));
            on_pyc.v(ti,1)=onsets.v(tempon(tind));
        end
    end

    sp_pyc=sp_py;
    for ti=2:length(sp_py.t)-1
        tempon=find(peaks.t>sp_py.t(ti-1) & peaks.t<sp_py.t(ti+1));
        if length(tempon)>0
            [~,tind]=min(abs(sp_py.t(ti)-peaks.t(tempon)));
            sp_pyc.t(ti,1)=peaks.t(tempon(tind));
            sp_pyc.v(ti,1)=peaks.v(tempon(tind));
        end
    end

    ext_on_ind=find(diff(on_pyc.t)==0);
    ext_sp_ind=find(diff(sp_pyc.t)==0);

    ext_all_ind=unique([ext_on_ind',ext_sp_ind']);

    on_pyc.t(ext_all_ind)=[];
    on_pyc.v(ext_all_ind)=[];
    on_pyc.i(ext_all_ind)=[];
    sp_pyc.t(ext_all_ind)=[];
    sp_pyc.v(ext_all_ind)=[];
    sp_pyc.i(ext_all_ind)=[];

    on_pym=on_py;
    sp_pym=sp_py;
    on_pym.t(ext_all_ind)=[];
    on_pym.v(ext_all_ind)=[];
    on_pym.i(ext_all_ind)=[];
    sp_pym.t(ext_all_ind)=[];
    sp_pym.v(ext_all_ind)=[];
    sp_pym.i(ext_all_ind)=[];

    new_on=eval(['on_',mode_on]);
    new_sp=eval(['sp_',mode_sp]);
end

function feat_data = calc_Xb200(peaks, onsets, fs, sig_data, up)

% eliminate any nans (which represent ectopics which have been removed)
peaks.t = peaks.t(~isnan(peaks.t));
peaks.v = peaks.v(~isnan(peaks.v));
onsets.t = onsets.t(~isnan(onsets.t));
onsets.v = onsets.v(~isnan(onsets.v));

% Find bw
bw.v = mean([onsets.v, peaks.v], 2);
bw.t = mean([onsets.t, peaks.t], 2);

% Find am
am.t = mean([onsets.t, peaks.t], 2);
am.v = [peaks.v - onsets.v];

% Get window modulation
% Define parameters
mu = 0;        % Mean of the Gaussian curve (center of the window)
sigma = 5;     % Standard deviation (controls the width of the curve)
dt=0.01;
win_for_amp=20;
t = -win_for_amp:dt:win_for_amp; % Time vector from 0 to 10 seconds with a step of 0.01 seconds

% Generate Gaussian curve using gaussmf
gaussian_curve = gaussmf(t, [sigma mu]);

prod=[];
for tmp_v=1:length(am.v)
    tmp_ind=find((am.t<am.t(tmp_v)+win_for_amp) &(am.t>am.t(tmp_v)-win_for_amp));
    tmp_amp=am.v(tmp_ind);
    tmp_sec=am.t(tmp_ind)-am.t(tmp_v)+win_for_amp;
    min_smpl=int32(min(tmp_sec)/dt)+1;
    max_smpl=int32(max(tmp_sec)/dt)+1;
    
    % Desired number of elements
    num_elem=(max_smpl-min_smpl)+1;
    
    % Equidistant time stamps
    equi_time = linspace(double(min_smpl*dt), double(max_smpl*dt), num_elem);
    
    % Interpolate values at equidistant time stamps using linear interpolation
    exp_values = interp1(tmp_sec, tmp_amp, equi_time, 'linear', 'extrap');

    tmp_gauss=gaussian_curve(min_smpl:max_smpl);

    prod(tmp_v)=sum(tmp_gauss.*exp_values)/double(num_elem);
end


% Normalise
feat_data.v = bw.v./nanmean(prod);
feat_data.t = bw.t;
feat_data.fs = fs;

end



function feat_data = calc_ageingindex(sp, on, fs, sig_data, up, bm)

% % eliminate any nans (which represent ectopics which have been removed)
% peaks.t = peaks.t(~isnan(peaks.t));
% peaks.v = peaks.v(~isnan(peaks.v));
% onsets.t = onsets.t(~isnan(onsets.t));
% onsets.v = onsets.v(~isnan(onsets.v));
% 
% % Find ageing index
% ageingindex.t = mean([onsets.t, peaks.t], 2);
% ageingindex.v = peaks.v;

% eliminate any nans (which represent ectopics which have been removed)
ind_AGI=AGI.indexs;

AGI = tmp_bm(~isnan(tmp_bm));
peaks.t = peaks(~isnan(peaks(ind_AGI)))/fs;
onsets.t = onsets(~isnan(onsets(ind_AGI)))/fs;


% Find ageing index
ageingindex.t = mean([onsets.t, peaks.t], 2);
ageingindex.v = AGI';

% Normalise
feat_data.v = ageingindex.v./nanmean(ageingindex.v);
feat_data.t = ageingindex.t;
feat_data.fs = fs;

end

function feat_data = calc_am(peaks, onsets, fs, sig_data, up)

% eliminate any nans (which represent ectopics which have been removed)
peaks.t = peaks.t(~isnan(peaks.t));
peaks.v = peaks.v(~isnan(peaks.v));
onsets.t = onsets.t(~isnan(onsets.t));
onsets.v = onsets.v(~isnan(onsets.v));

% Find am
am.t = mean([onsets.t, peaks.t], 2);
am.v = [peaks.v - onsets.v];

% Normalise
feat_data.v = am.v./nanmean(am.v);
feat_data.t = am.t;
feat_data.fs = fs;

end

function feat_data = calc_bw(peaks, onsets, fs, sig_data, up)

% eliminate any nans (which represent ectopics which have been removed)
peaks.t = peaks.t(~isnan(peaks.t));
peaks.v = peaks.v(~isnan(peaks.v));
onsets.t = onsets.t(~isnan(onsets.t));
onsets.v = onsets.v(~isnan(onsets.v));

% Find bw
bw.v = mean([onsets.v, peaks.v], 2);
bw.t = mean([onsets.t, peaks.t], 2);

% Find am
am.t = mean([onsets.t, peaks.t], 2);
am.v = [peaks.v - onsets.v];

% Normalise
feat_data.v = bw.v./nanmean(am.v);
feat_data.t = bw.t;
feat_data.fs = fs;

end


function feat_data = calc_bwm(peaks, onsets, fs, sig_data, up)

% eliminate any nans (which represent ectopics which have been removed)
peaks.t = peaks.t(~isnan(peaks.t));
peaks.v = peaks.v(~isnan(peaks.v));
onsets.t = onsets.t(~isnan(onsets.t));
onsets.v = onsets.v(~isnan(onsets.v));

% Find bwm
bwm.t = mean([onsets.t(2:end), onsets.t(1:(end-1))], 2);
bwm.v = nan(length(peaks.t)-1,1);
for s = 1 : (length(onsets.t)-1)
    rel_sig_els = sig_data.t >= onsets.t(s) & sig_data.t < onsets.t(s+1);
    bwm.v(s) = mean(sig_data.v(rel_sig_els));
end

% Find am
am.t = mean([onsets.t, peaks.t], 2);
am.v = [peaks.v - onsets.v];

% Normalise
feat_data.v = bwm.v./nanmean(am.v);
feat_data.t = bwm.t;
feat_data.fs = fs;

end

function feat_data = calc_fm(peaks, onsets, fs, sig_data, up)

% find fm
fm.v = [peaks.t(2:end) - peaks.t(1:(end-1))]/fs;
fm.t = mean([peaks.t(2:end), peaks.t(1:(end-1))], 2);

% eliminate any nans (which represent ectopics which have been removed)
fm.t = fm.t(~isnan(fm.t));
fm.v = fm.v(~isnan(fm.v));

% Normalise
feat_data.v = fm.v./nanmean(fm.v);
feat_data.t = fm.t;
feat_data.fs = fs;

end

function feat_data = calc_pk(peaks, onsets, fs, sig_data, up)

% eliminate any nans (which represent ectopics which have been removed)
peaks.t = peaks.t(~isnan(peaks.t));
peaks.v = peaks.v(~isnan(peaks.v));

feat_data = peaks;
feat_data.fs = fs;

end

function feat_data = calc_on(peaks, onsets, fs, sig_data, up)

% eliminate any nans (which represent ectopics which have been removed)
onsets.t = onsets.t(~isnan(onsets.t));
onsets.v = onsets.v(~isnan(onsets.v));

feat_data = onsets;
feat_data.fs = fs;

end

function feat_data = calc_qrsW(peaks, onsets, fs, sig_data, up)

% eliminate any nans (which represent ectopics which have been removed)
peaks.t = peaks.t(~isnan(peaks.t));
peaks.v = peaks.v(~isnan(peaks.v));

% Identify timings of QRS onsets and ends
[qrs_ends, qrs_onsets] = identify_qrs_timings(peaks, sig_data, up);

% QRS width
feat_data.t = mean([qrs_onsets.t(:), qrs_ends.t(:)], 2);
feat_data.v = qrs_ends.t - qrs_onsets.t;
feat_data.fs = fs;

end

function feat_data = calc_qrS(peaks, onsets, fs, sig_data, up)

% QRS slopes as described in:
% Lázaro, J. et al., 2014. Electrocardiogram Derived Respiratory Rate from QRS Slopes and R-Wave Angle. Annals of Biomedical Engineering, 42(10), pp.2072–83. Available at: http://www.ncbi.nlm.nih.gov/pubmed/25118665

% eliminate any nans (which represent ectopics which have been removed)
peaks.t = peaks.t(~isnan(peaks.t));
peaks.v = peaks.v(~isnan(peaks.v));
onsets.t = onsets.t(~isnan(onsets.t));
onsets.v = onsets.v(~isnan(onsets.v));

% Identify Q points
[~, qrs_onsets] = identify_qrs_timings(peaks, sig_data, up);

% Identify max upslopes between Q and R points
max_upslopes = identify_max_upslopes(peaks, qrs_onsets, sig_data, up);

% QR slopes
feat_data = max_upslopes;
feat_data.fs = fs;

end

function feat_data = calc_rsS(peaks, onsets, fs, sig_data, up)

% QRS slopes as described in:
% Lázaro, J. et al., 2014. Electrocardiogram Derived Respiratory Rate from QRS Slopes and R-Wave Angle. Annals of Biomedical Engineering, 42(10), pp.2072–83. Available at: http://www.ncbi.nlm.nih.gov/pubmed/25118665

% eliminate any nans (which represent ectopics which have been removed)
peaks.t = peaks.t(~isnan(peaks.t));
peaks.v = peaks.v(~isnan(peaks.v));
onsets.t = onsets.t(~isnan(onsets.t));
onsets.v = onsets.v(~isnan(onsets.v));

% Identify S points
[qrs_ends, ~] = identify_qrs_timings(peaks, sig_data, up);

% Identify max downslopes between Q and R points
max_downslopes = identify_max_downslopes(peaks, qrs_ends, sig_data, up);

% QR slopes
feat_data = max_downslopes;
feat_data.fs = fs;

end

function feat_data = calc_Rang(peaks, onsets, fs, sig_data, up)

% QRS angles as described in:
% Lázaro, J. et al., 2014. Electrocardiogram Derived Respiratory Rate from QRS Slopes and R-Wave Angle. Annals of Biomedical Engineering, 42(10), pp.2072–83. Available at: http://www.ncbi.nlm.nih.gov/pubmed/25118665

% eliminate any nans (which represent ectopics which have been removed)
peaks.t = peaks.t(~isnan(peaks.t));
peaks.v = peaks.v(~isnan(peaks.v));

% Identify timings of QRS onsets and ends
[qrs_ends, qrs_onsets] = identify_qrs_timings(peaks, sig_data, up);

% Identify max upslopes between Q and R points
max_upslopes = identify_max_upslopes(peaks, qrs_onsets, sig_data, up);

% Identify max downslopes between Q and R points
max_downslopes = identify_max_downslopes(peaks, qrs_ends, sig_data, up);

% identify angles between max upslopes and downslopes
if length(max_upslopes) ~= length(max_downslopes)
    error('Check this')
end

% Calculate QRS angles
qrs_angs = atan(abs( (max_upslopes.v - max_downslopes.v) ./ (0.4*(6.25 + (max_upslopes.v.*max_downslopes.v) ))));

% QRS angles
feat_data.t = mean([max_upslopes.t(:), max_downslopes.t(:)], 2);
feat_data.v = qrs_angs;
feat_data.fs = fs;

end

function feat_data = calc_qrsA(peaks, onsets, fs, sig_data, up)

sig_data.v = sig_data.v(:);
sig_data.t = sig_data.t(:);

% eliminate any nans (which represent ectopics which have been removed)
peaks.t = peaks.t(~isnan(peaks.t));
peaks.v = peaks.v(~isnan(peaks.v));

% Identify timings of QRS onsets and ends
[qrs_ends, qrs_onsets] = identify_qrs_timings(peaks, sig_data, up);

% QRS area
qrs_area.v = nan(length(qrs_onsets.t),1);
for s = 1 : length(qrs_onsets.t)
    rel_sig_els = qrs_onsets.i(s) : qrs_ends.i(s); rel_sig_els = rel_sig_els(:);
    % baseline
    baseline.t = sig_data.t(rel_sig_els);
    baseline.i = rel_sig_els;
    baseline.v = interp1([qrs_onsets.i(s), qrs_ends.i(s)], [sig_data.v(qrs_onsets.i(s)), sig_data.v(qrs_ends.i(s))], baseline.i);
    % baseline-corrected area
    qrs_area.v(s) = sum(sig_data.v(rel_sig_els)- baseline.v(:))/fs;
end

feat_data.t = mean([qrs_onsets.t(:), qrs_ends.t(:)], 2);
feat_data.v = qrs_area.v;
feat_data.fs = fs;

end

function [qrs_ends, qrs_onsets] = identify_qrs_timings(peaks, sig_data, up)

% Identify timings of QRS onsets and ends
thresh = up.paramSet.q_s_thresh;   % no of secs to search
[qrs_onsets.t, qrs_onsets.i, qrs_ends.t, qrs_ends.i] = deal(nan(length(peaks.t),1));
counter_no = 0;
for s = 1 : length(peaks.t)
    counter_no = counter_no+1;
    % Identify search regions
    rel_sig_els_on = sig_data.t < peaks.t(s) & sig_data.t > (peaks.t(s) - thresh);
    %[~,temp_on] = min(sig_data.v(rel_sig_els_on));
    rel_sig_els_end = sig_data.t > peaks.t(s) & sig_data.t < (peaks.t(s) + thresh);
    %[~,temp_end] = min(sig_data.v(rel_sig_els_end));
    
    if ~isempty(rel_sig_els_on) && ~isempty(rel_sig_els_end)
        % Find minimum closest to R-wave in each search region
        rel_data = sig_data.v(rel_sig_els_on);
        deriv = diff(rel_data);
        temp_troughs_on = deriv(2:end) >0 & deriv(1:(end-1)) <=0; temp_troughs_on = [0; temp_troughs_on(:); 0];
        rel_data = sig_data.v(rel_sig_els_end);
        deriv = diff(rel_data);
        temp_troughs_end = deriv(2:end) >0 & deriv(1:(end-1)) <=0; temp_troughs_end = [0; temp_troughs_end(:); 0];
        if sum(temp_troughs_on) ~= 0 & sum(temp_troughs_end) ~= 0
            qrs_onsets.i(counter_no) = max(find(temp_troughs_on))+ min(find(rel_sig_els_on))-1;
            qrs_onsets.t(counter_no) = sig_data.t(qrs_onsets.i(counter_no));
            qrs_ends.i(counter_no) = min(find(temp_troughs_end))+ min(find(rel_sig_els_end))-1;
            qrs_ends.t(counter_no) = sig_data.t(qrs_ends.i(counter_no));
        end
    end
    clear temp_troughs_on temp_troughs_end deriv rel_data rel_sig_els_on rel_sig_els_end
end

% eliminate remaining nans
for var = {'qrs_onsets', 'qrs_ends'}
    eval([var{1,1} '.i = ' var{1,1} '.i(~isnan(' var{1,1} '.i));']);
    eval([var{1,1} '.t = ' var{1,1} '.t(~isnan(' var{1,1} '.t));']);
end

% check that none of the QRS complexes overlap. If they do then delete the
% later one of the two:

temp = qrs_onsets.i(2:end) - qrs_ends.i(1:(end-1));
overlap_log = temp < 0;
if isempty(qrs_onsets.i)
    els_to_keep = [];
else
    els_to_keep = [1; find(~overlap_log(:)) + 1];
end

for var = {'qrs_onsets', 'qrs_ends'}
    eval([var{1,1} '.i = ' var{1,1} '.i(els_to_keep);']);
    eval([var{1,1} '.t = ' var{1,1} '.t(els_to_keep);']);
end

end

function max_upslopes = identify_max_upslopes(peaks, qrs_onsets, sig_data, up)

%% Identify peaks and QRS onsets

% Identify peaks for which a corresponding onset was found
good_peak_els = false(size(peaks.t));
for s = 1 : length(peaks.t)
    rel_onsets = qrs_onsets.t < peaks.t(s) & qrs_onsets.t >= peaks.t(s)-up.paramSet.q_s_thresh;
    if sum(rel_onsets)==1
        good_peak_els(s) = true;
    end
end

if sum(good_peak_els) ~= length(qrs_onsets.t)
    error('Check this')
end

% Eliminate peaks for which no corresponding end was found
peaks.t = peaks.t(good_peak_els);
peaks.v = peaks.v(good_peak_els);

% Find first derivative
der.v = sig_data.v(2:end) - sig_data.v(1:end-1);
der.t = sig_data.t(2:end)- (1/(2*sig_data.fs));

% Identify timings of upslopes
thresh = 0.004;  % no secs either side of max derivative to include when finding slope
no_els_either_side = thresh/(1/sig_data.fs); clear thresh
[max_upslopes.t, max_upslopes.v] = deal(nan(length(peaks.t),1));
for s = 1 : length(peaks.t)
    
    % Identify search region
    rel_der_els = sig_data.t < peaks.t(s) & sig_data.t > qrs_onsets.t(s);
    
    if ~isempty(rel_der_els)
        % Find time of maximum first derivative within search region
        rel_data = der.v(rel_der_els);
        [~, temp_el] = max(rel_data); clear rel_data
        found_rel_der_els = find(rel_der_els); 
        max_el = found_rel_der_els(temp_el); clear temp_el found_rel_der_els
        rel_els = (max_el - no_els_either_side) : (max_el + no_els_either_side);
        rel_data.y = sig_data.v(rel_els);
        rel_data.x = 1:length(rel_els); clear rel_els
        x_bar = mean(rel_data.x);
        y_bar = mean(rel_data.y);
        max_upslopes.t(s) = der.t(max_el); clear max_el
        max_upslopes.v(s) = sum((rel_data.x - x_bar).*(rel_data.y - y_bar))/sum((rel_data.x - x_bar).^2); clear y_bar x_bar rel_data
        
    end
    clear rel_der_els
end
clear s no_els_either_side

% eliminate remaining nans
for var = {'t', 'v'}
    eval(['max_upslopes.' var{1,1} ' = max_upslopes.' var{1,1} '(~isnan(max_upslopes.' var{1,1} '));']);
end

end

function max_downslopes = identify_max_downslopes(peaks, qrs_ends, sig_data, up)

%% Identify peaks and QRS ends

% Identify peaks for which a corresponding end was found
good_peak_els = false(size(peaks.t));
for s = 1 : length(peaks.t)
    rel_ends = qrs_ends.t > peaks.t(s) & qrs_ends.t <= peaks.t(s)+up.paramSet.q_s_thresh;
    if sum(rel_ends)==1
        good_peak_els(s) = true;
    end
end

if sum(good_peak_els) ~= length(qrs_ends.t)
    error('Check this')
end

% Eliminate peaks for which no corresponding end was found
peaks.t = peaks.t(good_peak_els);
peaks.v = peaks.v(good_peak_els);

%% Find downslopes

% Find first derivative
der.v = sig_data.v(2:end) - sig_data.v(1:end-1);
der.t = sig_data.t(2:end)- (1/(2*sig_data.fs));

% Identify downslopes
thresh = 0.004;  % no secs either side of max derivative to include when finding slope
no_els_either_side = thresh/(1/sig_data.fs); clear thresh
[max_downslopes.t, max_downslopes.v] = deal(nan(length(peaks.t),1));
for s = 1 : length(peaks.t)
    
    % Identify search region
    rel_der_els = sig_data.t > peaks.t(s) & sig_data.t < qrs_ends.t(s);
    
    if ~isempty(rel_der_els)
        % Find time of maximum first derivative within search region
        rel_data = der.v(rel_der_els);
        [~, temp_el] = min(rel_data); clear rel_data
        found_rel_der_els = find(rel_der_els); 
        min_el = found_rel_der_els(temp_el); clear temp_el found_rel_der_els
        rel_els = (min_el - no_els_either_side) : (min_el + no_els_either_side);
        rel_data.y = sig_data.v(rel_els);
        rel_data.x = 1:length(rel_els); clear rel_els
        x_bar = mean(rel_data.x);
        y_bar = mean(rel_data.y);
        max_downslopes.t(s) = der.t(min_el); clear min_el
        max_downslopes.v(s) = sum((rel_data.x - x_bar).*(rel_data.y - y_bar))/sum((rel_data.x - x_bar).^2); clear y_bar x_bar rel_data
        
    end
    clear rel_der_els
end
clear s no_els_either_side

% eliminate remaining nans
for var = {'t', 'v'}
    eval(['max_downslopes.' var{1,1} ' = max_downslopes.' var{1,1} '(~isnan(max_downslopes.' var{1,1} '));']);
end

end

function feat_data = calc_pca(peaks, onsets, fs, sig_data, up)
%PCA extracts a respiratory signal using kernel PCA as described in:
% Widjaja, D. et al., 2012. Application of kernel principal component analysis for single-lead-ECG-derived respiration. IEEE Transactions on Biomedical Engineering, 59(4), pp.1169–76. Available at: http://www.ncbi.nlm.nih.gov/pubmed/22438200.
%

%% This uses the LS-SVMlab toolbox
% The toolbox (v.1.8) can be downloaded from:
%    http://www.esat.kuleuven.be/sista/lssvmlab/
%
% The following is an excerpt from the website:
%
% The LS-SVMlab software is made available for non commercial research
% purposes only under the GNU General Public License. However,
% notwithstanding any provision of the GNU General Public License,
% LS-SVMlab software may not be used for commercial purposes without
% explicit written permission after contacting LS-SVMlab@esat.kuleuven.be.

% download the scripts if they're not in the search path
curr_dir = mfilename('fullpath'); curr_dir = curr_dir(1:end-3);
filepath = [curr_dir, 'LSSVMlabv1_8_R2009b_R2011a.zip'];
if ~exist(filepath, 'file')
    % download zip file
    url = 'http://www.esat.kuleuven.be/sista/lssvmlab/downloads/LSSVMlabv1_8_R2009b_R2011a.zip';
    downloadedfilename = websave(filepath,url);
    % unzip zip file
    unzip(downloadedfilename, curr_dir)
end
addpath(genpath(curr_dir))

%% Step 1: Input Matrix
% determine no of samples within either (i) 60 ms of QRS spike, or (ii) 300
% ms of ppg pulse:
if strcmp(sig_data.wave_type, 'ekg')
    beats = peaks;
    % eliminate any nans (which represent ectopics which have been removed)
    beats.t = beats.t(~isnan(beats.t));
    beats.v = beats.v(~isnan(beats.v));
    no_samps = round(0.06*fs);        % from: Widjaja, D. et al., 2012. Application of kernel principal component analysis for single-lead-ECG-derived respiration. IEEE Transactions on Biomedical Engineering, 59(4), pp.1169–76. Available at: http://www.ncbi.nlm.nih.gov/pubmed/22438200.
    int = 1;   % keep all samples
    tot_no_samps = (1+(no_samps*2));
else
    beats = onsets;                       % use onsets because they seem to be more consistent than peaks.
    % eliminate any nans (which represent ectopics which have been removed)
    temp = diff(beats.t);
    temp = temp(~isnan(temp));
    no_secs = 0.5*median(temp); clear temp  % consider the median beat length.
    beats.t = beats.t(~isnan(beats.t));
    beats.v = beats.v(~isnan(beats.v));
    no_samps = round(no_secs*fs);         % adapated from: Madhav, K.V. et al., 2010. Estimation of respiratory rate from principal components of photoplethysmographic signals. In 2010 IEEE EMBS Conference on Biomedical Engineering and Sciences (IECBES). IEEE, pp. 311–314. Available at: http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5742251 [Accessed November 12, 2013].
    int = round(fs/50);   % keep only those samples which are required to sample at 50 Hz (chosen because it's more than twice the highest frequency content retained (20 Hz))
    tot_no_samps = (1+(floor(2*no_samps/int)));
end

% construct input matrix of beats
input_mat = nan(length(beats.t), tot_no_samps);
for beat_no = 1 : length(beats.t)
    beat_el = find(sig_data.t == beats.t(beat_no));
    rel_els = (beat_el - no_samps) : int : (beat_el + no_samps);
    if rel_els(1) <= 0 || rel_els(end) > length(sig_data.v)
        continue
    end
    input_mat(beat_no,:) = sig_data.v(rel_els);
end
rel_t = beats.t(~isnan(input_mat(:,1)));
input_mat = input_mat(~isnan(input_mat(:,1)),:);

%% Step 2: Define omega squared values
omega_s_hat = size(input_mat, 2)*mean(var(input_mat));
omega_s = linspace((omega_s_hat/100), (omega_s_hat*100), 1000);

%% Step 3: Apply kPCA for each omega squared value
diffs = nan(length(omega_s),1);
[eigval, eigvec] = deal(cell(length(omega_s),1));
for omega_s_no = 1 : length(omega_s)
    [eigval{omega_s_no}, eigvec{omega_s_no}] = kpca(input_mat, 'RBF_kernel', omega_s(omega_s_no));
    %% Step 4: Select optimal omega squared
    % Here the difference between the first eigenvalue and the sum of the
    % remainder is maximised.
    first_eigval = eigval{omega_s_no}(1);
    sum_remaining_eigval = sum(eigval{omega_s_no}(2:end));
    diffs(omega_s_no) = first_eigval - sum_remaining_eigval;
    % break once you've passed the first maximum to save processing time.
    if omega_s_no>1 && diffs(omega_s_no) < diffs(omega_s_no-1)
        break
    end
end

[~, rel_el] = max(diffs);
rel_omega_s = omega_s(rel_el);
rel_eigvec = eigvec{rel_el};

%% Step 5: Reconstruction using the first eigenvector
first_eigvec = rel_eigvec(1,:);
Ximgn = preimage_rbf(input_mat,rel_omega_s,first_eigvec'); % ,B,type,npcs,maxIts

%% Step 6: Definition of respiratory signal
feat_data.v = Ximgn(:,1);
good_els = ~isnan(feat_data.v);
feat_data.v = detrend(feat_data.v(good_els));
feat_data.t = rel_t(good_els);
feat_data.fs = fs;

% %% Plot
% plot(sig_data.t, sig_data.v), hold on,
% plot(beats.t, beats.v, '.r'),
% plot(respWave.t, 3000*respWave.v)

end

function feat_data = calc_pulW(peaks, onsets, fs, sig_data, up)
%pulW extracts PPG pulse width using the methodology described in:
% Lázaro Plaza, J., 2015. Non-invasive techniques for respiratory information extraction based on pulse photoplethysmogram and electrocardiogram. Universidad Zaragoza.
%
nu = 0.05;

fs = sig_data.fs;
sig_data.v = sig_data.v(:)';
%% Step 1: LPF data
% LPF below 5 Hz
filt_characteristics.Fpass = 6.9;  % in HZ
filt_characteristics.Fstop = 4;  % in HZ   (18.36 and 23 provide a -3 dB cutoff of 20 Hz)
filt_characteristics.Dpass = 0.05;
filt_characteristics.Dstop = 0.01;
s_filt.v = elim_vhfs(sig_data.v, fs, filt_characteristics);
s_filt.t = sig_data.t;

%% Step 2: Find derivative of LPF data
s_der.t = [nan, s_filt.t(2:end)];
s_der.v = [nan, s_filt.v(2:end) - s_filt.v(1:(end-1))];
s_der2.v = [s_der.v(2:end) - s_der.v(1:(end-1)), nan];

[no, ne] = deal(nan(length(peaks.t),1));
for peak_no = 1 : length(peaks.t)
    if isnan(peaks.t(peak_no))
        continue
    end
    peak_el = find(s_der.t == peaks.t(peak_no));    % the peak in the original signal (before filtering)
    lower_el = peak_el - round(0.3*fs);             % 0.3s before the peak in the original signal
    if lower_el < 1
        continue
    end
    %% Step 3: Find nv
    % nv are the locations of the max derivatives before each pulse peak
    % (within 0.3 s before the peak)
    [~, temp] = max(s_der.v(lower_el:peak_el));
    nv = temp - 1 + lower_el;              % nv = the max upslope of the filtered signal during the 0.3 s before the peak
    clear temp peak_el
    
    %% Step 4: Identify search interval in which to look for onset
    int_els = lower_el:nv;                  % interval is from 0.3 s before the peak until the max upslope of the filtered signal.
    clear lower_el
    
    %% Step 5: See which criterion is met
    % evaluate C1
    thresh = nu*s_der.v(nv);
    if sum(s_der.v(int_els) <= thresh)
        c = 1;
    else
        % evaluate C2
        min_log = s_der.v(int_els(2:(end-1))) < s_der.v(int_els(1:(end-2))) & ...
            s_der.v(int_els(2:(end-1))) < s_der.v(int_els(3:end));
        if sum(min_log)
            c = 2;
        else
            c = 3;
        end
    end
    clear thresh
    
    %% Step 6: Find pulse onset according to criteria
    
    if c == 1
        [~, temp] = min(abs( s_der.v(int_els)- (nu*s_der.v(nv)) ));
        no(peak_no) = int_els(1) - 1 + temp;
    elseif c == 2
        no(peak_no) = int_els(1) + find(min_log, 1, 'last');
    else
        [~, temp] = min( s_der.v(int_els) );
        no(peak_no) = int_els(1) - 1 + temp;
    end
    clear min_log c temp int_els nv
    
    
    
    %% Repeat for pulse ends
    peak_el = find(s_der.t == peaks.t(peak_no));    % the peak in the original signal (before filtering)
    upper_el = peak_el + round(0.3*fs);             % 0.3s before the peak in the original signal
    if upper_el > length(s_filt.t)
        continue
    end
    %% Step 3: Find nd
    % nd are the locations of the min derivatives after each pulse peak
    % (within 0.3 s after the peak)
    [~, temp] = min(s_der.v(peak_el:upper_el));
    nd = temp - 1 + peak_el;              % nd = the min downslope of the filtered signal during the 0.3 s after the peak
    clear temp peak_el
    
    %% Step 4: Identify search interval in which to look for onset
    int_els = nd:upper_el;                  % interval is from the min downslope of the filtered signal until 0.3 s after the peak.
    clear upper_el
    
    %% Step 5: See which criterion is met
    % evaluate C1
    thresh = nu*s_der.v(nd);
    if sum(s_der.v(int_els) >= thresh)
        c = 1;
    else
        % evaluate C2
        min_log = s_der.v(int_els(2:(end-1))) > s_der.v(int_els(1:(end-2))) & ...
            s_der.v(int_els(2:(end-1))) > s_der.v(int_els(3:end));
        if sum(min_log)
            c = 2;
        else
            c = 3;
        end
    end
    clear thresh
    
    %% Step 6: Find pulse onset according to criteria
    
    if c == 1
        [~, temp] = min(abs( s_der.v(int_els)- (nu*s_der.v(nd)) ));
        ne(peak_no) = int_els(1) - 1 + temp;
    elseif c == 2
        ne(peak_no) = int_els(1) + find(min_log, 1, 'first');
    else
        [~, temp] = max( s_der.v(int_els) );
        ne(peak_no) = int_els(1) - 1 + temp;
    end
    clear min_log c temp int_els nd
    
end
clear peak_no


%     plot(sig_data.t, sig_data.v), hold on,
%     plot(s_filt.t, s_filt.v)
%     plot(s_filt.t, s_der2.v./range(s_der2.v))
%     plot(peaks.t, peaks.v, 'sk'),
%     plot(sig_data.t(no), sig_data.v(no), '*k'),
%     plot(sig_data.t(ne), sig_data.v(ne), 'ok')
%     close all

%% Step 7: Definition of respiratory signal
widths = (ne-no)./fs;
good_els  = ~isnan(widths);

feat_data.v = detrend(widths(good_els)); clear temp
feat_data.t = peaks.t(good_els); clear good_els
feat_data.fs = fs;


end