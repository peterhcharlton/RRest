function run_vortal_determinants_analysis
% RUN_VORTAL_DETERMINANTS_ANALYSIS performs statistical analyses of the 'vortal_factors'
% results.
%
% This will produce the results reported in the following publication:
%
%           Charlton P.H. et al. Extraction of respiratory signals from the 
%           electrocardiogram and photoplethysmogram: technical and physiological
%           determinants, Physiological Measurement, 38(5), 2017
%           DOI: https://doi.org/10.1088/1361-6579/aa670e
%
%	Inputs:
%		results         this script performs analysis of the outputs of
%                       RRest('vortal_factors'). Therefore, the results files
%                       from this operation must be available.
%
%	Outputs:
%       Results tables and figures, saved in:
%           ... \vortal_factors\Analysis_files\Publication_Results\
%       (note that these are numbered differently to the figures in the full text)
%           
%   Further Information:
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
up = setup_universal_params;

%% Load win_data file
win_data = load_win_data_file(up);

% %% Perform analysis of respiratory modulation indices
% resp_mod_inds(win_data, up);   % not used in final publication

%% Fig 0: Sampling Freqs
compare_samp_freqs(win_data, up);

%% Fig 1: Compare mods at different PPG sites
compare_ppg_site_diff_mods(win_data, up);

%% Fig 2: Compare lab and clinical mods
compare_lab_clinical_diff_mods(win_data, up);

%% Fig 3: Compare young and elderly
compare_young_and_elderly_box_mods(win_data, up);

%% Fig 4: Compare male and female
compare_gender_box_mods(win_data, up);
%compare_gender_box_mods_sub_groups(win_data, up);   % not used in final publication

%% Fig 5: Physiological factors
compare_rr_mods(win_data, up);
% compare_rr_mods_ind_subjs(win_data, up);   % not used in final publication
% compare_hr_mods(win_data, up);   % not used in final publication
% compare_hrrr_mods(win_data, up);   % not used in final publication

%% Fig 6: ECG vs PPG
compare_ecg_and_ppg_mods(win_data, up);

end

function up = setup_universal_params

fprintf('\n--- Creating Universal Parameters ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% FILE PATHS TO BE SPECIFIED %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% (requires editing) %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Please note that your system may require the slashes in file paths to be
% of the opposite direction. In which case, change the following:
up.paths.slash_direction = '\';     % usually a backslash for Windows, forward slash for Linux

% Specify path of data root folder. For instance, if the results of
% RRest('vortal_ext') are saved at:
%     'C:\Documents\Data\vortal_ext\Analysis_files\Results' ,
% then the data root folder should be specified as:
%     'C:\Documents\Data\'
up.paths.root_folder = 'C:\Documents\Data\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% LOAD PREVIOUS UNIVERSAL PARAMETERS %%%%%%%%
%%%%%%%%%%%%%%%% (no editing required) %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Provide licence details
provide_licence_details
% setup paths
up_filename = 'up';
dataset_name = 'vortal_factors';
up.paths.root_data_folder = [up.paths.root_folder, dataset_name, up.paths.slash_direction];
up.paths.data_save_folder = [up.paths.root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Component_Data', up.paths.slash_direction];
up_path = [up.paths.data_save_folder, up_filename];
% load up
load(up_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% ADD PARAMETERS SPECIFIC TO PUBLICATION %%%%%%%%
%%%%%%%%%%%%%%%% (no editing required) %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

up.pub.save_folder = [up.paths.root_data_folder, 'Analysis_files', up.paths.slash_direction, 'Publication_Results', up.paths.slash_direction];
if ~exist(up.pub.save_folder, 'dir')
    mkdir(up.pub.save_folder)
end

up.pub.rel_index = 'CCp';
up.pub.plot_lims = [-0.5, 0.5];
up.pub.y_ticks = -0.5:0.25:0.5;
up.pub.alpha = 0.05;
up.pub.equip_type = 'clin';  % clin or lab
close all  % incase any figs are open

end

function provide_licence_details

licence_details = ['\n\n run_vortal_determinants_analysis  Copyright (C) 2017  King''s College London',...
    '\n This program comes with ABSOLUTELY NO WARRANTY', ... 
    '\n and is available under the GNU public license.', ...
    '\n For details see the accompanying LICENSE.txt file.\n\n'];

fprintf(licence_details)

end

function win_data = load_win_data_file(up)

fprintf('\n--- Loading data');

% load win_data file
load([up.paths.data_save_folder, up.paths.filenames.win_data], 'win_data')

% eliminate non-sensical algorithms
% bad_els = (win_data.m_xb == 7 & ~win_data.ecg_log) | ...
%     (win_data.m_xb == 8 & ~win_data.ecg_log) | ...
%     win_data.m_xb == 10 | ...
%     win_data.m_xa == 4 | ...
%     ~win_data.comb_log | ...
%     win_data.subj == 10 | ...
%     win_data.ekgclinDS25Hz_log;

%bad_els = (~win_data.comb_log | win_data.ekgclinDS25Hz_log | win_data.m_xa == 4);
bad_els = (~win_data.comb_log | win_data.m_xa == 4);

fields = fieldnames(win_data);
for field_no = 1 : length(fields)
    eval(['win_data.' fields{field_no} ' = win_data.' fields{field_no} '(~bad_els);']);
end

%win_data = rmfield(win_data, 'ekgclinDS25Hz_log');

end

function resp_mod_inds(win_data, up)

% specify modulation measurement indices
inds = {'CCp', 'MSC', 'SNR'};
% specify RR estimation domains
rrest_doms = {'all'}; %, 'time', 'freq'};

% setup figure
lwidth = 2;
lbl_ftsize = 14;
paper_size = [200, 200, 1400, 600];
h_fig = figure('Position',paper_size);

% calculate errors
errors = abs(win_data.ref - win_data.est);
errors(~win_data.comb_log) = nan;

% identify errors which are relevant to each domain
for dom_no = 1 : length(rrest_doms)
    if strcmp(rrest_doms{dom_no}, 'time')
        dom_abr = 'Et';
        temp_rel_dom = ~cellfun(@isempty, strfind(win_data.alg_names, dom_abr));
        eval(['rel_dom.' rrest_doms{dom_no} ' = temp_rel_dom;']);
    elseif strcmp(rrest_doms{dom_no}, 'freq')
        dom_abr = 'Ef';
        temp_rel_dom = ~cellfun(@isempty, strfind(win_data.alg_names, dom_abr));
        eval(['rel_dom.' rrest_doms{dom_no} ' = temp_rel_dom;']);
    else
        rel_dom.all = true(length(win_data.subj),1);
    end
end

% cycle through indices
for ind_no = 1 : length(inds)
    curr_ind = inds{ind_no};
    eval(['rel_meas = win_data.' curr_ind ';']);
    
    % specify subplot
    subplot(1, length(inds), ind_no)
    
    % cycle through possible RR estimation domains
    for dom_no = 1 : length(rrest_doms)
        curr_dom = rrest_doms{dom_no};
        % eliminate elements corresponding to the other domain
        eval(['temp_rel_dom = rel_dom.' rrest_doms{dom_no} ';']);
        % identify errors corresponding to each bin of measure values
        bins.l = 0:39;
        bins.u = bins.l+1;
        [meas_vals.v, meas_vals.lq, meas_vals.uq, meas_n] = deal(nan(length(bins.l),1));
        for bin_no = 1 : length(bins.l)
            meas_n(bin_no) = sum(~isnan(rel_meas(errors>=bins.l(bin_no) & errors<bins.u(bin_no) & temp_rel_dom)));
            meas_vals.v(bin_no) = nanmedian(rel_meas(errors>=bins.l(bin_no) & errors<bins.u(bin_no) & temp_rel_dom));
            meas_vals.uq(bin_no) = quantile(rel_meas(errors>=bins.l(bin_no) & errors<bins.u(bin_no)), 0.25);
            meas_vals.lq(bin_no) = quantile(rel_meas(errors>=bins.l(bin_no) & errors<bins.u(bin_no)), 0.75);
        end
        % plot
        plot(meas_vals.v, 'LineWidth', lwidth), hold on,
        plot(meas_vals.lq, '--b', 'LineWidth', lwidth)
        plot(meas_vals.uq, '--b', 'LineWidth', lwidth)
        meas_n_stat = meas_n/max(meas_n);
        plot(meas_n_stat, 'r'),
        meas_n_stat(meas_n_stat>0.01) = nan;
        plot(meas_n_stat, '.k')
        hold on
        clear meas_vals
    end
    
    % add annotations to plot
    xlabel('Absolute Error [bpm]', 'FontSize', lbl_ftsize)
    ylabel(curr_ind, 'FontSize', lbl_ftsize)
    ylim([0 1])
    %legend(rrest_doms)
    
end
PrintFigs(gcf, paper_size(3:4)/100, [up.pub.save_folder, 'resp_mod_inds'])
end

function compare_samp_freqs(win_data, up)

fprintf('\n--- Fig 0');

% extract relevant index measures
eval(['rel_meas = win_data.' up.pub.rel_index ';']);

% specify cohorts
cohorts = {'all', 'elderly', 'young'};

% cycle through each cohort
for cohort_no = 1 : length(cohorts)
    curr_cohort = cohorts{cohort_no};
    % identify relevant rows for this cohort (taking only ppg rows)
    switch curr_cohort
        case 'all'
            rel_cohort_els = true(length(win_data.ecg_log),1);
        case 'young'
            rel_cohort_els = win_data.young_log;
        case 'elderly'
            rel_cohort_els = win_data.elderly_log;
    end
    
    % cycle through each sig
    for overall_sig = {'ekg', 'ppg'}
        
        % specify signal(s) of interest
        if strcmp(overall_sig, 'ppg')
            sig = 'ppgclin';
            control_sig = 'ppgclinDS125Hz';
        elseif strcmp(overall_sig, 'ekg')
            sig = 'ekgclin';
            control_sig = 'ekgclinDS500Hz';
        end
        
        % identify downsampled versions of signal
        temp = fieldnames(win_data);
        ds_sig_names = temp(~cellfun(@isempty, strfind(temp, sig)));  clear temp
        ds_sig_labels = ds_sig_names;
        for sig_no = 1 : length(ds_sig_labels)
            if isempty(strfind(ds_sig_labels{sig_no}, 'DS'))
                if strcmp(overall_sig, 'ppg')
                    ds_sig_labels{sig_no} = [ds_sig_labels{sig_no}(1:(end-4)), 'DS125Hz'];
                elseif strcmp(overall_sig, 'ekg')
                    ds_sig_labels{sig_no} = [ds_sig_labels{sig_no}(1:(end-4)), 'DS500Hz'];
                end
            else
                ds_sig_labels{sig_no} = ds_sig_labels{sig_no}(1:(end-4));
            end
        end
        
        % identify rr cutoffs
        if strcmp(overall_sig, 'ppg')
            order = [1,6,5,4,3,2,7];
            ds_sig_labels = ds_sig_labels(order);
            ds_sig_names = ds_sig_names(order);
        elseif strcmp(overall_sig, 'ekg')
            order = [1,5,4,3,2,7,6];
            ds_sig_labels = ds_sig_labels(order);
            ds_sig_names = ds_sig_names(order);
        end
        pops = ds_sig_labels;
        
        counter_no = 0; res = nan(1000,1000); [groups.mod, groups.sig, groups.sig] = deal(cell(0));
        for pop_no = 1:length(pops)
            
            % identify elements which are relevant to this signal
            eval(['rel_pop_els = win_data.' ds_sig_names{pop_no} ' & ~isnan(rel_meas);']);
            
            % specify modulations
            xb_mods = unique(win_data.m_xb(rel_pop_els)); xb_mods = xb_mods(~isnan(xb_mods));
            xa_mods = unique(win_data.m_xa(rel_pop_els)); xa_mods = xa_mods(~isnan(xa_mods));
            no_mods = length(xa_mods)+length(xb_mods);
            x_labs = cell(no_mods,1);
            for mod_no = 1 : no_mods
                if mod_no <= length(xa_mods)
                    x_labs{mod_no} = ['XA' num2str(xa_mods(mod_no))];
                else
                    x_labs{mod_no} = ['XB' num2str(xb_mods(mod_no-length(xa_mods)))];
                end
            end
            poss_mods = [xa_mods; xb_mods];
            
            for mod_no = 1 : length(poss_mods)
                curr_mod = poss_mods(mod_no);
                % identify relevant elements for this mod
                if mod_no <= length(xa_mods)
                    rel_mod_els = win_data.m_xa == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els & rel_cohort_els;
                else
                    rel_mod_els = win_data.m_xb == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els & rel_cohort_els;
                end
                % identify index measures for these elements
                temp_ind_meas = rel_meas(rel_mod_els);
                % identify win nos for these elements
                temp_win_no = win_data.win_no(rel_mod_els);
                % look at each subject
                temp_subj = win_data.subj(rel_mod_els);
                subjs = unique(temp_subj);
                subj_med_ind_meas = nan(1,length(subjs));
                for subj_no = 1:length(subjs)
                    curr_subj = subjs(subj_no);
                    subj_wins = unique(temp_win_no(temp_subj == curr_subj));
                    temp_meas = nan(length(subj_wins),1);
                    for win_no = 1 : length(subj_wins)
                        curr_win = subj_wins(win_no);
                        cand_ind_meas = unique(temp_ind_meas(temp_win_no == curr_win & temp_subj == curr_subj));
                        if length(cand_ind_meas) == 1
                            temp_meas(win_no,1) = cand_ind_meas;
                        else
                            error('Check this')
                        end
                    end
                    subj_med_ind_meas(subj_no) = median(temp_meas);
                end
                counter_no = counter_no+1;
                groups.sig{counter_no,1} = pops{pop_no};
                groups.mod{counter_no,1} = ['$' x_labs{mod_no}(1), '_{', x_labs{mod_no}(2:end), '}$'];   % for the latex interpreter
                res(counter_no,1:length(subj_med_ind_meas)) = subj_med_ind_meas;
                clear subj_med_ind_meas subj_no curr_subj subjs temp_subj
                clear curr_mod rel_mod_els
            end
            clear mod_no temp_ind_meas
            
        end
        res = res(1:length(groups.sig),1:length(unique(win_data.subj(rel_cohort_els)))); res = res';
        
        % adjust res for values at original sampling freq
        y_medians = nan(length(poss_mods),1);
        for mod_no = 1 : (length(poss_mods))
            if strcmp(overall_sig, 'ekg')
                young_mod_els = ~cellfun(@isempty, strfind(groups.sig, 'DS500')) & strcmp(groups.mod, groups.mod{mod_no});
            else
                young_mod_els = ~cellfun(@isempty, strfind(groups.sig, 'DS125')) & strcmp(groups.mod, groups.mod{mod_no});
            end
            temp_data = res(:,young_mod_els);
            y_medians(mod_no) = nanmedian(temp_data(:));
        end
        for col_no = 1 : size(res, 2)
            mod_no = find(strcmp(groups.mod, groups.mod{col_no}));
            res(:,col_no) = res(:,col_no) - y_medians(mod_no(1));
        end
        
        % statistical analysis
        [mod_stats.p, mod_stats.z] = deal(nan(size(res,2),1));
        for mod_no = 1 : size(res,2)
            control_el = strcmp(groups.sig, control_sig) & strcmp(groups.mod, groups.mod{mod_no});
            control_data = res(:,control_el);
            comparison_data = res(:,mod_no);
            rel_els = ~isnan(comparison_data) & ~isnan(control_data);
            control_data = control_data(rel_els);
            comparison_data = comparison_data(rel_els);
            [mod_stats.p(mod_no),~,stats] = signrank(comparison_data, control_data, 'method', 'exact');
            % check to see if stats.zval exists:
            if ~sum(strcmp(fieldnames(stats), 'zval')) && ~isequal(control_data, comparison_data)
                [~,~,stats2] = signrank(comparison_data, control_data, 'method', 'approximate');
                stats.zval = stats2.zval;
            else
                stats.zval = 0;
            end
            mod_stats.z(mod_no) = stats.zval;
        end
        
        % correction for multiple comparisons
        sig_diffs = correct_multiple_comparisons(mod_stats, up);
        
        % setup figure
        lwidth = 2;
        lbl_ftsize = 14;
        paper_size = [200, 200, 1000, 600];
        figure('Position',paper_size);
        
        % subplot 1
        subplot('Position', [0.22, 0.2, 0.75, 0.75])
        
        % make boxplot
        ylims = [-0.5, 0.5];
        new_labels = groups.mod(1:length(poss_mods));
        offset_temp = 3;
        new_labels = [repmat({''}, [length(unique(groups.sig))-offset_temp-1,length(new_labels)]); new_labels'; repmat({''}, [offset_temp,length(new_labels)])]; new_labels = new_labels(:);
        boxplot(res,{groups.mod,groups.sig},'colors','k','factorgap',[4 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
        hold on
        % add x-axis
        hline = refline([0 0]);
        hline.Color = 'k';
        % colour in boxes
        h = findobj(gca,'Tag','Box');
        [~, order] = sort(groups.mod);
        box_diffs = sig_diffs(order);
        box_ages = repmat(pops, [1, length(poss_mods)]);
        if strcmp(overall_sig, 'ppg')
            for j=1:length(h)
                if strcmp(box_ages(length(h)+1-j), ds_sig_labels{1})
                    prr1 = patch(get(h(j),'XData'),get(h(j),'YData'),0.9*[1,1,1]); % light grey
                elseif strcmp(box_ages(length(h)+1-j), ds_sig_labels{2})
                    prr2 = patch(get(h(j),'XData'),get(h(j),'YData'),0.8*[1,1,1]); % light grey
                elseif strcmp(box_ages(length(h)+1-j), ds_sig_labels{3})
                    prr3 = patch(get(h(j),'XData'),get(h(j),'YData'),0.7*[1,1,1]); % light grey
                elseif strcmp(box_ages(length(h)+1-j), ds_sig_labels{4})
                    prr4 = patch(get(h(j),'XData'),get(h(j),'YData'),0.6*[1,1,1]); % dark grey
                elseif strcmp(box_ages(length(h)+1-j), ds_sig_labels{5})
                    prr5 = patch(get(h(j),'XData'),get(h(j),'YData'),0.5*[1,1,1]); % dark grey
                elseif strcmp(box_ages(length(h)+1-j), ds_sig_labels{6})
                    prr6 = patch(get(h(j),'XData'),get(h(j),'YData'),0.4*[1,1,1]); % dark grey
                elseif strcmp(box_ages(length(h)+1-j), ds_sig_labels{7})
                    prr7 = patch(get(h(j),'XData'),get(h(j),'YData'),0.3*[1,1,1]); % dark grey
                end
            end
        else
            for j=1:length(h)
                if strcmp(box_ages(length(h)+1-j), ds_sig_labels{1})
                    prr1 = patch(get(h(j),'XData'),get(h(j),'YData'),0.9*[1,1,1]); % light grey
                elseif strcmp(box_ages(length(h)+1-j), ds_sig_labels{2})
                    prr2 = patch(get(h(j),'XData'),get(h(j),'YData'),0.8*[1,1,1]); % light grey
                elseif strcmp(box_ages(length(h)+1-j), ds_sig_labels{3})
                    prr3 = patch(get(h(j),'XData'),get(h(j),'YData'),0.7*[1,1,1]); % light grey
                elseif strcmp(box_ages(length(h)+1-j), ds_sig_labels{4})
                    prr4 = patch(get(h(j),'XData'),get(h(j),'YData'),0.6*[1,1,1]); % dark grey
                elseif strcmp(box_ages(length(h)+1-j), ds_sig_labels{5})
                    prr5 = patch(get(h(j),'XData'),get(h(j),'YData'),0.5*[1,1,1]); % dark grey
                elseif strcmp(box_ages(length(h)+1-j), ds_sig_labels{6})
                    prr6 = patch(get(h(j),'XData'),get(h(j),'YData'),0.4*[1,1,1]); % dark grey
                elseif strcmp(box_ages(length(h)+1-j), ds_sig_labels{7})
                    prr7 = patch(get(h(j),'XData'),get(h(j),'YData'),0.3*[1,1,1]); % dark grey
                end
            end
        end
        % replot boxplot, highlighting significant differences
        boxplot(res,{groups.mod,groups.sig},'colors','k','factorgap',[4 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
        for j=1:length(h)
            if box_diffs(length(h)+1-j) == 1
                lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
                xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
                ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
                rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',[0.5,0.5,1], 'LineWidth',lwidth)
            elseif box_diffs(length(h)+1-j) == -1
                lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
                xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
                ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
                rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',0.8*[1,1,0], 'LineWidth',lwidth)
            end
        end
        % colour in median line
        h = findobj(gca,'tag','Median');
        set(h,'Color','r')
        % y-axis label
        ylab = ylabel({'CC', 'relative to', 'median CC', 'at highest', 'sampling', 'frequency', }, 'FontSize', lbl_ftsize, 'Rotation', 0);
        set(ylab, 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0], 'VerticalAlignment', 'middle');
        set(gca, 'FontSize', lbl_ftsize, 'TickLabelInterpreter', 'latex')
        box off
        grid on
        set(gca,'XGrid', 'off')
        ylim(up.pub.plot_lims), set(gca,'YTick', up.pub.y_ticks)
        
        % create legend
        leg_labels = ds_sig_labels;
        for s = 1 : length(leg_labels)
            leg_labels{s} = leg_labels{s}(10:end);
            if length(leg_labels{s}) == 3
                leg_labels{s} = ['    ', leg_labels{s}];
            elseif length(leg_labels{s}) == 4
                leg_labels{s} = ['  ', leg_labels{s}];
            end
        end
        if strcmp(overall_sig, 'ppg')
            leg_h = legend([prr1, prr2, prr3, prr4, prr5, prr6, prr7], leg_labels{1}, leg_labels{2}, leg_labels{3}, leg_labels{4}, leg_labels{5}, leg_labels{6}, leg_labels{7});
        else
            leg_h = legend([prr1, prr2, prr3, prr4, prr5, prr6, prr7], leg_labels{1}, leg_labels{2}, leg_labels{3}, leg_labels{4}, leg_labels{5}, leg_labels{6}, leg_labels{7});
        end
        set(leg_h, 'FontSize', lbl_ftsize - 4, 'Position',[0.01,0.74,0.16,0.17])
        
        % Annotate highest freqs at which quality is significantly reduced:
        mods = unique(groups.mod);
        no_freqs = length(box_diffs)/length(mods);
        new_new_labels.mod_no = [];
        new_new_labels.str = cell(0);
        for mod_no = 1 : length(mods)
            curr_mod = strrep(mods{mod_no}, ' ', '');
            rel_box_diffs = box_diffs(1+no_freqs*(mod_no-1) : no_freqs*mod_no) ;
            if sum(rel_box_diffs) == 0
                continue
            else
                first_rel_box_diff = find(rel_box_diffs, 1, 'first');
                first_freq = strrep(leg_labels{first_rel_box_diff}, ' ', '');
                new_new_labels.mod_no(end+1) = mod_no;
                new_new_labels.str{end+1} = first_freq;
            end
        end
        left_hand_x_pos = 0.195; x_width = 0.77;
        annotation('textbox',[0 0.05 left_hand_x_pos .1],'String',{'First', 'reduction', 'in CC:'},'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
        x_width_per_mod = x_width/length(poss_mods);
        start_x_pos = left_hand_x_pos - 0.002;
        x_labs2 = cell(1,length(poss_mods));
        for mod_no = 1 : length(poss_mods)
            rel_box_diffs = box_diffs(1+no_freqs*(mod_no-1) : no_freqs*mod_no) ;
            if sum(new_new_labels.mod_no == mod_no) == 0
                x_labs2{mod_no} = '';    % NS
            else
                first_rel_box_diff = find(rel_box_diffs, 1, 'first');
                x_labs2{mod_no} = strrep(leg_labels{first_rel_box_diff}, ' ', '');
            end
            dim = [start_x_pos+(x_width_per_mod*(mod_no-1)) 0 .1 .1];
            str = x_labs2{mod_no};
            annotation('textbox',dim,'String',str,'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
        end
        
        % save plot
        filename = ['fig0_' overall_sig{1,1} '_' cohorts{cohort_no}];
        PrintFigs(gcf, paper_size(3:4)/100, [up.pub.save_folder, filename])
        
    end
    
end

end

function PrintFigs(h, paper_size, savepath)
set(gcf,'color','w');
set(h,'PaperUnits','inches');
set(h,'PaperSize', [paper_size(1), paper_size(2)]);
set(h,'PaperPosition',[0 0 paper_size(1) paper_size(2)]);
print(h,'-dpdf',savepath)
print(h,'-dpng',savepath)

% if you want .eps illustrations, then do as follows:
up.eps_figs = 1;
if up.eps_figs
    % you need to download 'export_fig' from:
    % http://uk.mathworks.com/matlabcentral/fileexchange/23629-export-fig
    export_fig_dir_path = 'C:\Documents\Google Drive\Work\Projects\PhD\Github\phd\Tools\Other Scripts\export_fig\altmany-export_fig-76bd7fa\';
    addpath(genpath(export_fig_dir_path))
    export_fig(savepath, '-eps')
end
close all
end

function compare_ppg_site_diff_mods(win_data, up)

fprintf('\n--- Fig 1');

% specify and extract modulation measurement
eval(['rel_meas = win_data.' up.pub.rel_index ';']);

% specify two signals of interest
sigs = {'ppgfraw', 'ppgeraw'};

% specify cohorts
cohorts = {'all', 'young', 'elderly'};

% cycle through each cohort
for cohort_no = 1 : length(cohorts)
    curr_cohort = cohorts{cohort_no};
    % identify relevant rows for this cohort (taking only ppg rows)
    switch curr_cohort
        case 'all'
            rel_cohort_els = ~win_data.ecg_log;
        case 'young'
            rel_cohort_els = win_data.young_log & ~win_data.ecg_log;
        case 'elderly'
            rel_cohort_els = win_data.elderly_log & ~win_data.ecg_log;
    end
    
    % specify subjects
    subjs = unique(win_data.subj(rel_cohort_els));
    
    % specify modulations
    xb_mods = unique(win_data.m_xb(rel_cohort_els)); xb_mods = xb_mods(~isnan(xb_mods));
    xa_mods = unique(win_data.m_xa(rel_cohort_els)); xa_mods = xa_mods(~isnan(xa_mods));
    no_mods = length(xa_mods)+length(xb_mods);
    x_labs = cell(no_mods,1);
    for mod_no = 1 : no_mods
        if mod_no <= length(xa_mods)
            x_labs{mod_no} = ['Xa' num2str(xa_mods(mod_no))];
        else
            x_labs{mod_no} = ['Xb' num2str(xb_mods(mod_no-length(xa_mods)))];
        end        
    end
    poss_mods = [xa_mods; xb_mods];
    
    % cycle through each modulation
    mod_meas = nan(length(poss_mods),length(subjs));
    for mod_no = 1 : length(poss_mods)
        curr_mod = poss_mods(mod_no);
        
        % extract relevant data for this cohort and modulation
        if mod_no <= length(xa_mods)
            rel_els = rel_cohort_els & win_data.m_xa == curr_mod;
        else
            rel_els = rel_cohort_els & win_data.m_xb == curr_mod;
        end
        temp_data.meas = rel_meas(rel_els);
        temp_data.subj = win_data.subj(rel_els);
        eval(['temp_data.sig1_log = win_data.' sigs{1} '_log(rel_els);']);
        eval(['temp_data.sig2_log = win_data.' sigs{2} '_log(rel_els);']);
        temp_data.win_no = win_data.win_no(rel_els);
        
        % cycle through each subj
        for subj_no = 1 : length(subjs)
            curr_subj = subjs(subj_no);
            subj_els = temp_data.subj == curr_subj;
            sig1_wins = unique(temp_data.win_no(subj_els & temp_data.sig1_log));
            sig2_wins = unique(temp_data.win_no(subj_els & temp_data.sig2_log));
            subj_wins = intersect(sig1_wins, sig2_wins);
            temp_meas = nan(length(subj_wins),2);
            for win_no = 1 : length(subj_wins)
                curr_win_no = subj_wins(win_no);
                sig1_win_els = subj_els & temp_data.sig1_log & temp_data.win_no == curr_win_no & ~isnan(temp_data.meas);
                sig2_win_els = subj_els & temp_data.sig2_log & temp_data.win_no == curr_win_no & ~isnan(temp_data.meas);
                if length(unique(temp_data.meas(sig1_win_els))) == 1
                    temp_meas(win_no,1) = unique(temp_data.meas(sig1_win_els));
                else
                    error('Check this')
                end
                if length(unique(temp_data.meas(sig2_win_els))) == 1
                    temp_meas(win_no,2) = unique(temp_data.meas(sig2_win_els));
                else
                    error('Check this')
                end
            end
            mod_meas(mod_no,subj_no) = median(temp_meas(:,1) - temp_meas(:,2));
        end
    end
    
    % statistical analysis
    [mod_stats.p, mod_stats.z] = deal(nan(length(poss_mods),1));
    for mod_no = 1 : length(poss_mods)
        rel_data = mod_meas(mod_no,:);
        [mod_stats.p(mod_no),~,stats] = signrank(rel_data, zeros(size(rel_data)), 'method', 'exact');
        % check to see if stats.zval exists:
        if ~sum(strcmp(fieldnames(stats), 'zval'))
            [~,~,stats2] = signrank(rel_data, zeros(size(rel_data)), 'method', 'approximate');
            stats.zval = stats2.zval;
        end
        mod_stats.z(mod_no) = stats.zval;
    end
    % correction for multiple comparisons
    sig_diffs = correct_multiple_comparisons(mod_stats, up);
    
    % Results table
    [results.n] = deal(nan(length(poss_mods),1));
    [results.med, results.lq, results.uq, results.mod, results.stat_sig] = deal(cell(length(poss_mods),1));
    results.mod = strrep(x_labs, ' ', '');
    results.stat_sig(sig_diffs == 1) = {'Fin'};
    results.stat_sig(sig_diffs == -1) = {'Ear'};
    for mod_no = 1 : length(poss_mods)
        rel_data = mod_meas(mod_no,:);
        results.n(mod_no) = sum(~isnan(rel_data));
        results.med{mod_no} = sprintf('%.2f', nanmedian(rel_data));
        results.lq{mod_no} = sprintf('%.2f', quantile(rel_data(~isnan(rel_data)), 0.25));
        results.uq{mod_no} = sprintf('%.2f', quantile(rel_data(~isnan(rel_data)), 0.75));
    end
    results.tbl = table(results.mod, results.n, results.stat_sig, results.med, results.lq, results.uq);
    eval(['results.' cohorts{cohort_no} ' = results.tbl;'])
    
    % setup figure
    lwidth = 4;
    lbl_ftsize = 14;
    paper_size = [200, 200, 1000, 500];
    figure('Position',paper_size);
    
    % subplot 1
    subplot('Position', [0.2, 0.2, 0.75, 0.75])
    
    % boxplot
    boxplot(mod_meas'), hold on
    % add x-axis
    hline = refline([0 0]);
    hline.Color = 'k';
    % colour in boxes
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        % colour in box
        patch(get(h(j),'XData'),get(h(j),'YData'),0.8*[1,1,1]);
        % coloured outline for statistically significant diffs
        if sig_diffs(length(h)+1-j) == 1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',[0.5,0.5,1], 'LineWidth',lwidth)
        elseif sig_diffs(length(h)+1-j) == -1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',0.8*[1,1,0], 'LineWidth',lwidth)
        end
    end
    % colour in median line
    h = findobj(gca,'tag','Median');
    set(h,'Color','r')
    % y-axis label
    ylab = ylabel({'median', '(CCfin -', 'CCear)'}, 'FontSize', lbl_ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
    set(gca, 'FontSize', lbl_ftsize)
    boxplot(mod_meas')
    % color in edges of boxes
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        set(h(j), 'Color', [0,0,0])
    end
    % add mod labels
    set(gca,'XTickLabel',x_labs, 'FontSize', lbl_ftsize)
    box off
    grid on
    set(gca,'XGrid', 'off')
    ylim(up.pub.plot_lims), set(gca,'YTick', up.pub.y_ticks)
    
    % subplot 2
    left_hand_x_pos = 0.18; x_width = 0.75;
    annotation('textbox',[0 0.05 left_hand_x_pos .1],'String','Significantly greater CC:','LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    x_width_per_mod = x_width/length(poss_mods);
    start_x_pos = left_hand_x_pos - 0.002;
    x_labs2 = cell(1,length(poss_mods));
    for s = 1 : length(poss_mods)
        if sig_diffs(s) == 0
            x_labs2{s} = '';   % NS
        elseif sig_diffs(s) == 1
            x_labs2{s} = 'finger';
        elseif sig_diffs(s) == -1
            x_labs2{s} = 'ear';
        end
        dim = [start_x_pos+(x_width_per_mod*(s-1)) 0 .1 .1];
        str = x_labs2{s};
        annotation('textbox',dim,'String',str,'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    end
    
    % save plot
    filename = ['fig1_' cohorts{cohort_no}];
    PrintFigs(gcf, paper_size(3:4)/100, [up.pub.save_folder, filename])
    
    eval(['tbl.' cohorts{cohort_no} ' = results.' cohorts{cohort_no} ';']);
    
end

% save results table
filename = 'tbl_1';
save([up.pub.save_folder, filename], 'tbl')

end

function compare_lab_clinical_diff_mods(win_data, up)

fprintf('\n--- Fig 2');

% specify and extract modulation measurement
eval(['rel_meas = win_data.' up.pub.rel_index ';']);

for sig_cat = {'ekg', 'ppg'}
    
    % specify two signals of interest
    switch sig_cat{1,1}
        case 'ppg'
            sigs = {'ppgclin', 'ppgfraw'};
        case 'ekg'
            sigs = {'ekgclin', 'ekgraw'};
    end
    
    % specify cohorts
    cohorts = {'all', 'young', 'elderly'};
    
    % cycle through each cohort
    for cohort_no = 1 : length(cohorts)
        curr_cohort = cohorts{cohort_no};
        switch curr_cohort
            case 'all'
                rel_cohort_els = true(length(win_data.subj),1);
            case 'young'
                rel_cohort_els = win_data.young_log;
            case 'elderly'
                rel_cohort_els = win_data.elderly_log;
        end
        
        % specify subjects
        subjs = unique(win_data.subj(rel_cohort_els));
        
        % specify modulations
        switch sig_cat{1,1}
            case 'ppg'
                xb_mods = unique(win_data.m_xb(rel_cohort_els & ~win_data.ecg_log)); xb_mods = xb_mods(~isnan(xb_mods));
                xa_mods = unique(win_data.m_xa(rel_cohort_els & ~win_data.ecg_log)); xa_mods = xa_mods(~isnan(xa_mods));
            case 'ekg'
                xb_mods = unique(win_data.m_xb(rel_cohort_els & win_data.ecg_log)); xb_mods = xb_mods(~isnan(xb_mods));
                xa_mods = unique(win_data.m_xa(rel_cohort_els & win_data.ecg_log)); xa_mods = xa_mods(~isnan(xa_mods));
        end
        no_mods = length(xa_mods)+length(xb_mods);
        x_labs = cell(no_mods,1);
        for mod_no = 1 : no_mods
            if mod_no <= length(xa_mods)
                x_labs{mod_no} = ['Xa' num2str(xa_mods(mod_no))];
            else
                x_labs{mod_no} = ['Xb' num2str(xb_mods(mod_no-length(xa_mods)))];
            end
        end
        poss_mods = [xa_mods; xb_mods];
        
        % cycle through each modulation
        mod_meas = nan(length(poss_mods),length(subjs));
        for mod_no = 1 : no_mods
            curr_mod = poss_mods(mod_no);
            
            % extract relevant data for this cohort and modulation
            if mod_no <= length(xa_mods)
                rel_els = rel_cohort_els & win_data.m_xa == curr_mod;
            else
                rel_els = rel_cohort_els & win_data.m_xb == curr_mod;
            end
            temp_data.meas = rel_meas(rel_els);
            temp_data.subj = win_data.subj(rel_els);
            eval(['temp_data.sig1_log = win_data.' sigs{1} '_log(rel_els);']);
            eval(['temp_data.sig2_log = win_data.' sigs{2} '_log(rel_els);']);
            temp_data.win_no = win_data.win_no(rel_els);
            
            % cycle through each subj
            for subj_no = 1 : length(subjs)
                curr_subj = subjs(subj_no);
                subj_els = temp_data.subj == curr_subj;
                sig1_wins = unique(temp_data.win_no(subj_els & temp_data.sig1_log));
                sig2_wins = unique(temp_data.win_no(subj_els & temp_data.sig2_log));
                subj_wins = intersect(sig1_wins, sig2_wins);
                temp_meas = nan(length(subj_wins),2);
                for win_no = 1 : length(subj_wins)
                    curr_win_no = subj_wins(win_no);
                    sig1_win_els = subj_els & temp_data.sig1_log & temp_data.win_no == curr_win_no & ~isnan(temp_data.meas);
                    sig2_win_els = subj_els & temp_data.sig2_log & temp_data.win_no == curr_win_no & ~isnan(temp_data.meas);
                    if length(unique(temp_data.meas(sig1_win_els))) == 1
                        temp_meas(win_no,1) = unique(temp_data.meas(sig1_win_els));
                    else
                        error('Check this')
                    end
                    if length(unique(temp_data.meas(sig2_win_els))) == 1
                        temp_meas(win_no,2) = unique(temp_data.meas(sig2_win_els));
                    else
                        error('Check this')
                    end
                end
                mod_meas(mod_no,subj_no) = median(temp_meas(:,1) - temp_meas(:,2));
            end
        end
        
        % statistical analysis
        [mod_stats.p, mod_stats.z] = deal(nan(length(poss_mods),1));
        for mod_no = 1 : length(poss_mods)
            rel_data = mod_meas(mod_no,:);
            [mod_stats.p(mod_no),~,stats] = signrank(rel_data, zeros(size(rel_data)), 'method', 'exact');
            % check to see if stats.zval exists:
            if ~sum(strcmp(fieldnames(stats), 'zval'))
                [~,~,stats2] = signrank(rel_data, zeros(size(rel_data)), 'method', 'approximate');
                stats.zval = stats2.zval;
            end
            mod_stats.z(mod_no) = stats.zval;
        end
        % correction for multiple comparisons
        sig_diffs = correct_multiple_comparisons(mod_stats, up);
        
        % Results table
        [results.n] = deal(nan(length(poss_mods),1));
        [results.med, results.lq, results.uq, results.mod, results.stat_sig] = deal(cell(length(poss_mods),1));
        results.mod = strrep(x_labs, ' ', '');
        results.stat_sig(sig_diffs == 1) = {'clin'};
        results.stat_sig(sig_diffs == -1) = {'lab'};
        for mod_no = 1 : length(poss_mods)
            rel_data = mod_meas(mod_no,:);
            results.n(mod_no) = sum(~isnan(rel_data));
            results.med{mod_no} = sprintf('%.2f', nanmedian(rel_data));
            results.lq{mod_no} = sprintf('%.2f', quantile(rel_data(~isnan(rel_data)), 0.25));
            results.uq{mod_no} = sprintf('%.2f', quantile(rel_data(~isnan(rel_data)), 0.75));
        end
        results.tbl = table(results.mod, results.n, results.stat_sig, results.med, results.lq, results.uq);
        eval(['results.' sig_cat{1,1} '.' cohorts{cohort_no} ' = results.tbl;'])
        
        
        % setup figure
        lwidth = 4;
        lbl_ftsize = 14;
        paper_size = [200, 200, 1000, 500];
        figure('Position',paper_size);
        
        % subplot 1
        subplot('Position', [0.2, 0.2, 0.75, 0.75])
        
        % boxplot
        boxplot(mod_meas'), hold on
        % add x-axis
        hline = refline([0 0]);
        hline.Color = 'k';
        % colour in boxes
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            % colour in box
            patch(get(h(j),'XData'),get(h(j),'YData'),0.8*[1,1,1]);
            % coloured outline for statistically significant diffs
            if sig_diffs(length(h)+1-j) == 1
                lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
                xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
                ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
                rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',[0.5,0.5,1], 'LineWidth',lwidth)
            elseif sig_diffs(length(h)+1-j) == -1
                lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
                xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
                ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
                rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',0.8*[1,1,0], 'LineWidth',lwidth)
            end
        end
        % colour in median line
        h = findobj(gca,'tag','Median');
        set(h,'Color','r')
        % y-axis label
        ylab = ylabel({'median', '(CCclin -', 'CClab)'}, 'FontSize', lbl_ftsize, 'Rotation', 0);
        set(ylab, 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
        set(gca, 'FontSize', lbl_ftsize)
        boxplot(mod_meas');
        % color in edges of boxes
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            set(h(j), 'Color', [0,0,0])
        end
        % add mod labels
        set(gca,'XTickLabel',x_labs, 'FontSize', lbl_ftsize)
        box off
        grid on
        set(gca,'XGrid', 'off')
        ylim(up.pub.plot_lims), set(gca,'YTick', up.pub.y_ticks)
        
        % subplot 2
        left_hand_x_pos = 0.18; x_width = 0.75;
        annotation('textbox',[0 0.05 left_hand_x_pos .1],'String','Significantly greater CC:','LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
        x_width_per_mod = x_width/length(poss_mods);
        start_x_pos = left_hand_x_pos - 0.002;
        x_labs2 = cell(1,length(poss_mods));
        for s = 1 : length(poss_mods)
            if sig_diffs(s) == 0
                x_labs2{s} = '';    % NS
            elseif sig_diffs(s) == 1
                x_labs2{s} = 'clin';
            elseif sig_diffs(s) == -1
                x_labs2{s} = 'lab';
            end
            dim = [start_x_pos+(x_width_per_mod*(s-1)) 0 .1 .1];
            str = x_labs2{s};
            annotation('textbox',dim,'String',str,'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
        end
        
        % save plot
        filename = ['fig2_' sig_cat{1,1} '_' cohorts{cohort_no}];
        PrintFigs(gcf, paper_size(3:4)/100, [up.pub.save_folder, filename])
        
        % store results table
        eval(['tbl.' sig_cat{1,1} '_' cohorts{cohort_no} ' = results.' sig_cat{1,1} '.' cohorts{cohort_no} ';']);
        
    end
    
end

% save results table
filename = 'tbl_2';
save([up.pub.save_folder, filename], 'tbl')

end

function compare_young_and_elderly_box_mods(win_data, up)

fprintf('\n--- Fig 3');

% extract relevant index measures
eval(['rel_meas = win_data.' up.pub.rel_index ';']);

% cycle through each sig
for overall_sig = {'ppg', 'ekg'}
    
    % specify signal(s) of interest
    if strcmp(up.pub.equip_type, 'clin')
        if strcmp(overall_sig, 'ppg')
            sig = 'ppgclin';
        elseif strcmp(overall_sig, 'ekg')
            sig = 'ekgclin';
        end
    elseif strcmp(up.pub.equip_type, 'lab')
        if strcmp(overall_sig, 'ppg')
            sig = 'ppgfraw';
        elseif strcmp(overall_sig, 'ekg')
            sig = 'ekgraw';
        end
    end
        
    
    % identify elements which are relevant to this signal
    eval(['rel_sig_els = win_data.' sig '_log;']);
    
    pops = {'young', 'elderly'};
    counter_no = 0; res = nan(1000,1000); [groups.mod, groups.age] = deal(cell(0)); 
    for pop_no = 1:length(pops)
        
        switch pops{pop_no}
            case 'all'
                rel_pop_els = true(length(win_data),1);
            case 'young'
                rel_pop_els = win_data.young_log;
            case 'elderly'
                rel_pop_els = win_data.elderly_log;
            case 'young_male'
                rel_pop_els = win_data.young_log & win_data.male_log;
            case 'young_female'
                rel_pop_els = win_data.young_log & win_data.female_log;
            case 'elderly_male'
                rel_pop_els = win_data.elderly_log & win_data.male_log;
            case 'elderly_female'
                rel_pop_els = win_data.elderly_log & win_data.female_log;
        end
        
        % specify modulations
        xb_mods = unique(win_data.m_xb(rel_sig_els)); xb_mods = xb_mods(~isnan(xb_mods));
        xa_mods = unique(win_data.m_xa(rel_sig_els)); xa_mods = xa_mods(~isnan(xa_mods));
        no_mods = length(xa_mods)+length(xb_mods);
        x_labs = cell(no_mods,1);
        for mod_no = 1 : no_mods
            if mod_no <= length(xa_mods)
                x_labs{mod_no} = ['Xa' num2str(xa_mods(mod_no))];
            else
                x_labs{mod_no} = ['Xb' num2str(xb_mods(mod_no-length(xa_mods)))];
            end
        end
        poss_mods = [xa_mods; xb_mods];
        
        for mod_no = 1 : length(poss_mods)
            curr_mod = poss_mods(mod_no);
            % identify relevant elements for this mod            
            if mod_no <= length(xa_mods)
                rel_mod_els = rel_sig_els & win_data.m_xa == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            else
                rel_mod_els = rel_sig_els & win_data.m_xb == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            end
            % identify index measures for these elements
            temp_ind_meas = rel_meas(rel_mod_els);
            % identify win nos for these elements
            temp_win_no = win_data.win_no(rel_mod_els);
            % look at each subject
            temp_subj = win_data.subj(rel_mod_els);
            subjs = unique(temp_subj);
            subj_med_ind_meas = nan(1,length(subjs));
            for subj_no = 1:length(subjs)
                curr_subj = subjs(subj_no);
                subj_wins = unique(temp_win_no(temp_subj == curr_subj));
                temp_meas = nan(length(subj_wins),1);
                for win_no = 1 : length(subj_wins)
                    curr_win = subj_wins(win_no);
                    cand_ind_meas = unique(temp_ind_meas(temp_win_no == curr_win & temp_subj == curr_subj));
                    if length(cand_ind_meas) == 1
                        temp_meas(win_no,1) = cand_ind_meas;
                    else
                        error('Check this')
                    end
                end
                subj_med_ind_meas(subj_no) = median(temp_meas);
            end
            counter_no = counter_no+1;
            groups.age{counter_no,1} = pops{pop_no};
            groups.mod{counter_no,1} = ['    ' x_labs{mod_no}];
            res(counter_no,1:length(subj_med_ind_meas)) = subj_med_ind_meas;
            clear subj_med_ind_meas subj_no curr_subj subjs temp_subj curr_mod rel_mod_els
        end
        clear mod_no temp_ind_meas
        
    end
    res = res(1:length(groups.age),:); res = res';
    
    % Results table
    [results.y.med, results.y.lq, results.y.uq, results.y.n, results.mod, results.stat_sig] = deal(cell(length(poss_mods),1));
    results.e = results.y;
    results.mod = strrep(x_labs, ' ', '');
    for mod_no = 1 : length(poss_mods)
        rel_data.y = res(:,strcmp(groups.mod, groups.mod{mod_no}) & strcmp(groups.age, 'young')); rel_data.y = rel_data.y(~isnan(rel_data.y));
        rel_data.e = res(:,strcmp(groups.mod, groups.mod{mod_no}) & strcmp(groups.age, 'elderly')); rel_data.e = rel_data.e(~isnan(rel_data.e));
        results.y.n{mod_no} = sprintf('%.2f', length(rel_data.y));
        results.e.n{mod_no} = sprintf('%.2f', length(rel_data.e));
        results.y.med{mod_no} = sprintf('%.2f', median(rel_data.y));
        results.y.lq{mod_no} = sprintf('%.2f', quantile(rel_data.y, 0.25));
        results.y.uq{mod_no} = sprintf('%.2f', quantile(rel_data.y, 0.75));
        results.e.med{mod_no} = sprintf('%.2f', median(rel_data.e));
        results.e.lq{mod_no} = sprintf('%.2f', quantile(rel_data.e, 0.25));
        results.e.uq{mod_no} = sprintf('%.2f', quantile(rel_data.e, 0.75));
    end
    
    % adjust res for young median values
    y_medians = nan(length(poss_mods),1);
    for mod_no = 1 : (length(poss_mods))
        young_mod_els = strcmp(groups.age, 'young') & strcmp(groups.mod, groups.mod{mod_no});
        y_medians(mod_no) = nanmedian(res(:,young_mod_els));        
    end
    for col_no = 1 : size(res, 2)
        mod_no = find(strcmp(groups.mod, groups.mod{col_no}));
        res(:,col_no) = res(:,col_no) - y_medians(mod_no(1));
    end
    
    % find any statistically significant differences
    [mod_stats.p, mod_stats.z] = deal(nan(length(poss_mods),1));
    for mod_no = 1 : length(poss_mods)
        rel_data.y = res(:,strcmp(groups.mod, groups.mod{mod_no}) & strcmp(groups.age, 'young')); rel_data.y = rel_data.y(~isnan(rel_data.y));
        rel_data.e = res(:,strcmp(groups.mod, groups.mod{mod_no}) & strcmp(groups.age, 'elderly')); rel_data.e = rel_data.e(~isnan(rel_data.e));
        [mod_stats.p(mod_no),~,stats] = ranksum(rel_data.y, rel_data.e, 'method', 'exact');
        % check to see if stats.zval exists:
        if ~sum(strcmp(fieldnames(stats), 'zval'))
            [~,~,stats2] = ranksum(rel_data.y, rel_data.e, 'method', 'approximate');
            stats.zval = stats2.zval;
        end
        mod_stats.z(mod_no) = stats.zval;
    end
    % correction for multiple comparisons
    sig_diffs = correct_multiple_comparisons(mod_stats, up);
    
    % significant differences
    results.stat_sig(sig_diffs == 1) = {'Young'};
    results.stat_sig(sig_diffs == -1) = {'Elderly'};
    results.tbl = table(results.mod, results.y.n, results.y.med, results.y.lq, results.y.uq, results.e.n, results.e.med, results.e.lq, results.e.uq, results.stat_sig);
    eval(['results.' overall_sig{1,1} ' = results.tbl;'])
        
    % setup figure
    lwidth = 3;
    lbl_ftsize = 14;
    paper_size = [200, 200, 1000, 600];
    figure('Position',paper_size);
    
    % subplot 1
    subplot('Position', [0.22, 0.2, 0.75, 0.75])
    
    % make boxplot
    ylims = [-0.5, 0.5];
    new_labels = groups.mod(1:length(poss_mods));
    new_labels = [new_labels'; repmat({''}, [1,length(new_labels)])]; new_labels = new_labels(:);
    boxplot(res,{groups.mod,groups.age},'colors','k','factorgap',[8 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
    hold on
    % add x-axis
    hline = refline([0 0]);
    hline.Color = 'k';
    % colour in boxes
    h = findobj(gca,'Tag','Box');
    box_diffs = [sig_diffs(:)'; -1*sig_diffs(:)']; box_diffs = box_diffs(:);
    box_ages = repmat({'young', 'elderly'}, [1, length(poss_mods)]);
    for j=1:length(h)
        if strcmp(box_ages(length(h)+1-j), 'young')
            py = patch(get(h(j),'XData'),get(h(j),'YData'),0.8*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'elderly')
            pe = patch(get(h(j),'XData'),get(h(j),'YData'),0.5*[1,1,1]); % dark grey
        end
    end
    % replot boxplot
    boxplot(res,{groups.mod,groups.age},'colors','k','factorgap',[8 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
    for j=1:length(h)
        if box_diffs(length(h)+1-j) == 1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',[0.5,0.5,1], 'LineWidth',lwidth)
        elseif box_diffs(length(h)+1-j) == -1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',0.8*[1,1,0], 'LineWidth',lwidth)
        end
    end
    % colour in median line
    h = findobj(gca,'tag','Median');
    set(h,'Color','r')
    % y-axis label
    ylab = ylabel({'CC', 'relative', 'to median', 'CCyoung'}, 'FontSize', lbl_ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0], 'VerticalAlignment', 'middle');
    set(gca, 'FontSize', lbl_ftsize)
    box off
    grid on
    set(gca,'XGrid', 'off')
    ylim(up.pub.plot_lims), set(gca,'YTick', up.pub.y_ticks)
    
    % subplot 2
    left_hand_x_pos = 0.195; x_width = 0.77;
    annotation('textbox',[0 0.025 left_hand_x_pos .1],'String',{'Significantly', 'larger CC:'},'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    x_width_per_mod = x_width/length(poss_mods);
    start_x_pos = left_hand_x_pos - 0.002;
    x_labs2 = cell(1,length(poss_mods));
    for s = 1 : length(poss_mods)
        if sig_diffs(s) == 0
            x_labs2{s} = '';    % NS
        elseif sig_diffs(s) == 1
            x_labs2{s} = 'young';
        elseif sig_diffs(s) == -1
            x_labs2{s} = 'elderly';
        end
        dim = [start_x_pos+(x_width_per_mod*(s-1)) 0 .1 .1];
        str = x_labs2{s};
        annotation('textbox',dim,'String',str,'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    end
    
    legend([py, pe], 'young', 'elderly', 'Location', 'best')
    
    % save plot
    filename = ['fig3_' overall_sig{1,1}];
    PrintFigs(gcf, paper_size(3:4)/100, [up.pub.save_folder, filename])
    
    % store results table
    eval(['tbl.' overall_sig{1,1} ' = results.' overall_sig{1,1} ';']);
    
end

% save results table
filename = 'tbl_3';
save([up.pub.save_folder, filename], 'tbl')

end

function compare_gender_box_mods(win_data, up)

fprintf('\n--- Fig 4');

% extract relevant index measures
eval(['rel_meas = win_data.' up.pub.rel_index ';']);

% cycle through each sig
for overall_sig = {'ppg', 'ekg'}
    
    % specify signal(s) of interest
    if strcmp(up.pub.equip_type, 'clin')
        if strcmp(overall_sig, 'ppg')
            sig = 'ppgclin';
        elseif strcmp(overall_sig, 'ekg')
            sig = 'ekgclin';
        end
    elseif strcmp(up.pub.equip_type, 'lab')
        if strcmp(overall_sig, 'ppg')
            sig = 'ppgfraw';
        elseif strcmp(overall_sig, 'ekg')
            sig = 'ekgraw';
        end
    end
    
    % identify elements which are relevant to this signal
    eval(['rel_sig_els = win_data.' sig '_log;']);
    
    pops = {'male', 'female'};
    counter_no = 0; res = nan(1000,1000); [groups.mod, groups.age] = deal(cell(0)); 
    for pop_no = 1:length(pops)
        
        switch pops{pop_no}
            case 'all'
                rel_pop_els = true(length(win_data),1);
            case 'male'
                rel_pop_els = win_data.male_log;
            case 'female'
                rel_pop_els = win_data.female_log;
        end
        
        % specify modulations
        xb_mods = unique(win_data.m_xb(rel_sig_els)); xb_mods = xb_mods(~isnan(xb_mods));
        xa_mods = unique(win_data.m_xa(rel_sig_els)); xa_mods = xa_mods(~isnan(xa_mods));
        no_mods = length(xa_mods)+length(xb_mods);
        x_labs = cell(no_mods,1);
        for mod_no = 1 : no_mods
            if mod_no <= length(xa_mods)
                x_labs{mod_no} = ['Xa' num2str(xa_mods(mod_no))];
            else
                x_labs{mod_no} = ['Xb' num2str(xb_mods(mod_no-length(xa_mods)))];
            end
        end
        poss_mods = [xa_mods; xb_mods];
        
        for mod_no = 1 : length(poss_mods)
            curr_mod = poss_mods(mod_no);
            % identify relevant elements for this mod            
            if mod_no <= length(xa_mods)
                rel_mod_els = rel_sig_els & win_data.m_xa == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            else
                rel_mod_els = rel_sig_els & win_data.m_xb == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            end
            % identify index measures for these elements
            temp_ind_meas = rel_meas(rel_mod_els);
            % identify win nos for these elements
            temp_win_no = win_data.win_no(rel_mod_els);
            % look at each subject
            temp_subj = win_data.subj(rel_mod_els);
            subjs = unique(temp_subj);
            subj_med_ind_meas = nan(1,length(subjs));
            for subj_no = 1:length(subjs)
                curr_subj = subjs(subj_no);
                subj_wins = unique(temp_win_no(temp_subj == curr_subj));
                temp_meas = nan(length(subj_wins),1);
                for win_no = 1 : length(subj_wins)
                    curr_win = subj_wins(win_no);
                    cand_ind_meas = unique(temp_ind_meas(temp_win_no == curr_win & temp_subj == curr_subj));
                    if length(cand_ind_meas) == 1
                        temp_meas(win_no,1) = cand_ind_meas;
                    else
                        error('Check this')
                    end
                end
                subj_med_ind_meas(subj_no) = median(temp_meas);
            end
            counter_no = counter_no+1;
            groups.age{counter_no,1} = pops{pop_no};
            groups.mod{counter_no,1} = ['    ' x_labs{mod_no}];
            res(counter_no,1:length(subj_med_ind_meas)) = subj_med_ind_meas;
            clear subj_med_ind_meas subj_no curr_subj subjs temp_subj
            clear curr_mod rel_mod_els
        end
        clear mod_no temp_ind_meas
        
    end
    res = res(1:length(groups.age),:); res = res';
    
    % adjust res for young median values
    m_medians = nan(length(poss_mods),1);
    for mod_no = 1 : (length(poss_mods))
        male_mod_els = strcmp(groups.age, 'male') & strcmp(groups.mod, groups.mod{mod_no});
        temp_data = res(:,male_mod_els);
        m_medians(mod_no) = nanmedian(temp_data(:));        
    end
    for col_no = 1 : size(res, 2)
        mod_no = find(strcmp(groups.mod, groups.mod{col_no}));
        res(:,col_no) = res(:,col_no) - m_medians(mod_no(1));
    end
    
    % find any statistically significant differences
    [mod_stats.p, mod_stats.z] = deal(nan(length(poss_mods),1));
    for mod_no = 1 : length(poss_mods)
        rel_data.m = res(:,strcmp(groups.mod, groups.mod{mod_no}) & strcmp(groups.age, 'male')); rel_data.m = rel_data.m(~isnan(rel_data.m));
        rel_data.f = res(:,strcmp(groups.mod, groups.mod{mod_no}) & strcmp(groups.age, 'female')); rel_data.f = rel_data.f(~isnan(rel_data.f));
        [mod_stats.p(mod_no),~,stats] = ranksum(rel_data.m, rel_data.f, 'method', 'exact');
        % check to see if stats.zval exists:
        if ~sum(strcmp(fieldnames(stats), 'zval'))
            [~,~,stats2] = ranksum(rel_data.m, rel_data.f, 'method', 'approximate');
            stats.zval = stats2.zval;
        end
        mod_stats.z(mod_no) = stats.zval;
    end
    % correction for multiple comparisons
    sig_diffs = correct_multiple_comparisons(mod_stats, up);
    sig_diffs_m = sig_diffs;
    sig_diffs_f = -1*sig_diffs;
    
    % setup figure
    lwidth = 3;
    lbl_ftsize = 14;
    paper_size = [200, 200, 1000, 600];
    figure('Position',paper_size);
    
    % subplot 1
    subplot('Position', [0.22, 0.2, 0.75, 0.75])
    
    % make boxplot
    ylims = [-0.5, 0.5];
    new_labels = groups.mod(1:length(poss_mods));
    new_labels = [new_labels'; repmat({''}, [1,length(new_labels)])]; new_labels = new_labels(:);
    boxplot(res,{groups.mod,groups.age},'colors','k','factorgap',[8 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
    hold on
    % add x-axis
    hline = refline([0 0]);
    hline.Color = 'k';
    % colour in boxes
    h = findobj(gca,'Tag','Box');
    box_diffs = [sig_diffs_m(:)'; sig_diffs_f(:)']; box_diffs = box_diffs(:);
    box_ages = repmat({'male', 'female'}, [1, length(poss_mods)]);
    for j=1:length(h)
        if strcmp(box_ages(length(h)+1-j), 'male')
            pm = patch(get(h(j),'XData'),get(h(j),'YData'),0.8*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'female')
            pf = patch(get(h(j),'XData'),get(h(j),'YData'),0.4*[1,1,1]); % light grey
        end
    end
    % replot boxplot
    boxplot(res,{groups.mod,groups.age},'colors','k','factorgap',[8 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
    for j=1:length(h)
        if box_diffs(length(h)+1-j) == 1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',[0.5,0.5,1], 'LineWidth',lwidth)
        elseif box_diffs(length(h)+1-j) == -1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',0.8*[1,1,0], 'LineWidth',lwidth)
        end
    end
    % colour in median line
    h = findobj(gca,'tag','Median');
    set(h,'Color','r')
    % y-axis label
    ylab = ylabel({'CC', 'relative', 'to median', 'CCmale'}, 'FontSize', lbl_ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0], 'VerticalAlignment', 'middle');
    set(gca, 'FontSize', lbl_ftsize)
    box off
    grid on
    set(gca,'XGrid', 'off')
    ylim(up.pub.plot_lims), set(gca,'YTick', up.pub.y_ticks)
    
    % subplot 2
    left_hand_x_pos = 0.195; x_width = 0.77;
    annotation('textbox',[0 0.025 left_hand_x_pos .1],'String',{'Significantly', 'larger CC:'},'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    x_width_per_mod = x_width/length(poss_mods);
    start_x_pos = left_hand_x_pos - 0.002;
    x_labs2 = cell(1,length(poss_mods));
    for s = 1 : length(poss_mods)
        if sig_diffs_f(s) == 0
            x_labs2{s} = '';    % NS
        elseif sig_diffs_f(s) == 1
            x_labs2{s} = 'female';
        elseif sig_diffs_f(s) == -1
            x_labs2{s} = 'male';
        end
        dim = [start_x_pos+(x_width_per_mod*(s-1)) 0 .1 .1];
        str = x_labs2{s};
        annotation('textbox',dim,'String',str,'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    end
    
    legend([pm, pf], 'male', 'female', 'Location', 'best')
    
    % save plot
    filename = ['fig4_' overall_sig{1,1}];
    PrintFigs(gcf, paper_size(3:4)/100, [up.pub.save_folder, filename])
    
end

end

function compare_gender_box_mods_sub_groups(win_data, up)

fprintf('\n--- Fig 4 subgroup');

% extract relevant index measures
eval(['rel_meas = win_data.' up.pub.rel_index ';']);

% cycle through each sig
for overall_sig = {'ekg', 'ppg'}
    
    % specify signal(s) of interest
    if strcmp(overall_sig, 'ppg')
        sig = 'ppgclin';
    elseif strcmp(overall_sig, 'ekg')
        sig = 'ekgclin';
    end
    
    % identify elements which are relevant to this signal
    eval(['rel_sig_els = win_data.' sig '_log;']);
    
    pops = {'young_male', 'young_female', 'elderly_male', 'elderly_female'};
    counter_no = 0; res = nan(1000,1000); [groups.mod, groups.age] = deal(cell(0)); 
    for pop_no = 1:length(pops)
        
        switch pops{pop_no}
            case 'all'
                rel_pop_els = true(length(win_data),1);
            case 'young'
                rel_pop_els = win_data.young_log;
            case 'elderly'
                rel_pop_els = win_data.elderly_log;
            case 'young_male'
                rel_pop_els = win_data.young_log & win_data.male_log;
            case 'young_female'
                rel_pop_els = win_data.young_log & win_data.female_log;
            case 'elderly_male'
                rel_pop_els = win_data.elderly_log & win_data.male_log;
            case 'elderly_female'
                rel_pop_els = win_data.elderly_log & win_data.female_log;
        end
        
        % specify modulations
        xb_mods = unique(win_data.m_xb(rel_sig_els)); xb_mods = xb_mods(~isnan(xb_mods));
        xa_mods = unique(win_data.m_xa(rel_sig_els)); xa_mods = xa_mods(~isnan(xa_mods));
        no_mods = length(xa_mods)+length(xb_mods);
        x_labs = cell(no_mods,1);
        for mod_no = 1 : no_mods
            if mod_no <= length(xa_mods)
                x_labs{mod_no} = ['Xa' num2str(xa_mods(mod_no))];
            else
                x_labs{mod_no} = ['Xb' num2str(xb_mods(mod_no-length(xa_mods)))];
            end
        end
        poss_mods = [xa_mods; xb_mods];
        
        for mod_no = 1 : length(poss_mods)
            curr_mod = poss_mods(mod_no);
            % identify relevant elements for this mod            
            if mod_no <= length(xa_mods)
                rel_mod_els = rel_sig_els & win_data.m_xa == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            else
                rel_mod_els = rel_sig_els & win_data.m_xb == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            end
            % identify index measures for these elements
            temp_ind_meas = rel_meas(rel_mod_els);
            % identify win nos for these elements
            temp_win_no = win_data.win_no(rel_mod_els);
            % look at each subject
            temp_subj = win_data.subj(rel_mod_els);
            subjs = unique(temp_subj);
            subj_med_ind_meas = nan(1,length(subjs));
            for subj_no = 1:length(subjs)
                curr_subj = subjs(subj_no);
                subj_wins = unique(temp_win_no(temp_subj == curr_subj));
                temp_meas = nan(length(subj_wins),1);
                for win_no = 1 : length(subj_wins)
                    curr_win = subj_wins(win_no);
                    cand_ind_meas = unique(temp_ind_meas(temp_win_no == curr_win & temp_subj == curr_subj));
                    if length(cand_ind_meas) == 1
                        temp_meas(win_no,1) = cand_ind_meas;
                    else
                        error('Check this')
                    end
                end
                subj_med_ind_meas(subj_no) = median(temp_meas);
            end
            counter_no = counter_no+1;
            groups.age{counter_no,1} = pops{pop_no};
            groups.mod{counter_no,1} = ['    ' x_labs{mod_no}];
            res(counter_no,1:length(subj_med_ind_meas)) = subj_med_ind_meas;
            clear subj_med_ind_meas subj_no curr_subj subjs temp_subj
            clear curr_mod rel_mod_els
        end
        clear mod_no temp_ind_meas
        
    end
    res = res(1:length(groups.age),:); res = res';
    
    % adjust res for young median values
    y_medians = nan(length(poss_mods),1);
    for mod_no = 1 : (length(poss_mods))
        young_mod_els = ~cellfun(@isempty, strfind(groups.age, 'young')) & strcmp(groups.mod, groups.mod{mod_no});
        temp_data = res(:,young_mod_els);
        y_medians(mod_no) = nanmedian(temp_data(:));        
    end
    for col_no = 1 : size(res, 2)
        mod_no = find(strcmp(groups.mod, groups.mod{col_no}));
        res(:,col_no) = res(:,col_no) - y_medians(mod_no(1));
    end
    
    % find any statistically significant differences
    [mod_stats_y.p, mod_stats_y.z, mod_stats_e.p, mod_stats_e.z] = deal(nan(length(poss_mods),1));
    for mod_no = 1 : length(poss_mods)
        rel_data.ym = res(:,strcmp(groups.mod, groups.mod{mod_no}) & strcmp(groups.age, 'young_male')); rel_data.ym = rel_data.ym(~isnan(rel_data.ym));
        rel_data.yf = res(:,strcmp(groups.mod, groups.mod{mod_no}) & strcmp(groups.age, 'young_female')); rel_data.yf = rel_data.yf(~isnan(rel_data.yf));
        rel_data.em = res(:,strcmp(groups.mod, groups.mod{mod_no}) & strcmp(groups.age, 'elderly_male')); rel_data.em = rel_data.em(~isnan(rel_data.em));
        rel_data.ef = res(:,strcmp(groups.mod, groups.mod{mod_no}) & strcmp(groups.age, 'elderly_female')); rel_data.ef = rel_data.ef(~isnan(rel_data.ef));
        [mod_stats_y.p(mod_no),~,stats] = ranksum(rel_data.ym, rel_data.yf, 'method', 'exact');
        % check to see if stats.zval exists:
        if ~sum(strcmp(fieldnames(stats), 'zval'))
            [~,~,stats2] = ranksum(rel_data.ym, rel_data.yf, 'method', 'approximate');
            stats.zval = stats2.zval;
        end
        mod_stats_y.z(mod_no) = stats.zval;
        [mod_stats_e.p(mod_no),~,stats] = ranksum(rel_data.em, rel_data.ef, 'method', 'exact');
        % check to see if stats.zval exists:
        if ~sum(strcmp(fieldnames(stats), 'zval'))
            [~,~,stats2] = ranksum(rel_data.em, rel_data.ef, 'method', 'approximate');
            stats.zval = stats2.zval;
        end
        mod_stats_e.z(mod_no) = stats.zval;
    end
    % correction for multiple comparisons
    sig_diffs_y = correct_multiple_comparisons(mod_stats_y, up);
    sig_diffs_e = correct_multiple_comparisons(mod_stats_e, up);
    
    % setup figure
    lwidth = 3;
    lbl_ftsize = 14;
    paper_size = [200, 200, 1000, 600];
    figure('Position',paper_size);
    
    % subplot 1
    subplot('Position', [0.22, 0.2, 0.75, 0.75])
    
    % make boxplot
    ylims = [-0.5, 0.5];
    new_labels = groups.mod(1:length(poss_mods));
    new_labels = [new_labels'; repmat({''}, [3,length(new_labels)])]; new_labels = new_labels(:);
    boxplot(res,{groups.mod,groups.age},'colors','k','factorgap',[8 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
    hold on
    % add x-axis
    hline = refline([0 0]);
    hline.Color = 'k';
    % colour in boxes
    h = findobj(gca,'Tag','Box');
    box_diffs = [sig_diffs_y(:)'; -1*sig_diffs_y(:)'; sig_diffs_e(:)'; -1*sig_diffs_e(:)']; box_diffs = box_diffs(:);
    box_ages = repmat({'young_male', 'young_female', 'elderly_male', 'elderly_female'}, [1, length(poss_mods)]);
    for j=1:length(h)
        if strcmp(box_ages(length(h)+1-j), 'young_male')
            pym = patch(get(h(j),'XData'),get(h(j),'YData'),0.9*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'young_female')
            pyf = patch(get(h(j),'XData'),get(h(j),'YData'),0.7*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'elderly_male')
            pem = patch(get(h(j),'XData'),get(h(j),'YData'),0.5*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'elderly_female')
            pef = patch(get(h(j),'XData'),get(h(j),'YData'),0.3*[1,1,1]); % dark grey
        end
    end
    % replot boxplot
    boxplot(res,{groups.mod,groups.age},'colors','k','factorgap',[8 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
    for j=1:length(h)
        if box_diffs(length(h)+1-j) == 1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',[0.5,0.5,1], 'LineWidth',lwidth)
        elseif box_diffs(length(h)+1-j) == -1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',0.8*[1,1,0], 'LineWidth',lwidth)
        end
    end
    % colour in median line
    h = findobj(gca,'tag','Median');
    set(h,'Color','r')
    % y-axis label
    ylab = ylabel({'CC', 'relative', 'to median', 'CCyoung'}, 'FontSize', lbl_ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0], 'VerticalAlignment', 'middle');
    set(gca, 'FontSize', lbl_ftsize)
    box off
    grid on
    set(gca,'XGrid', 'off')
    ylim(up.pub.plot_lims), set(gca,'YTick', up.pub.y_ticks)
    
    % subplot 2
    left_hand_x_pos = 0.195; x_width = 0.77;
    annotation('textbox',[0 0.025 left_hand_x_pos .1],'String',{'Significantly', 'larger CC:'},'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    x_width_per_mod = x_width/length(poss_mods);
    start_x_pos = left_hand_x_pos - 0.002;
    x_labs2 = cell(1,length(poss_mods));
    for s = 1 : length(poss_mods)
        if sig_diffs_y(s) == 0 && sig_diffs_e(s) == 0
            x_labs2{s} = '';    % NS
        elseif sig_diffs_y(s) == 1 && sig_diffs_e(s) == 0
            x_labs2{s} = {'young', 'male'};
        elseif sig_diffs_y(s) == -1 && sig_diffs_e(s) == 0
            x_labs2{s} = {'young', 'female'};
        elseif sig_diffs_e(s) == 1 && sig_diffs_y(s) == 0
            x_labs2{s} = {'elderly', 'male'};
        elseif sig_diffs_e(s) == -1 && sig_diffs_y(s) == 0
            x_labs2{s} = {'elderly', 'female'};
        elseif sig_diffs_e(s) == -1 && sig_diffs_y(s) == -1
            x_labs2{s} = {'both', 'female'};
        elseif sig_diffs_e(s) == 1 && sig_diffs_y(s) == 1
            x_labs2{s} = {'both', 'male'};
        elseif sig_diffs_e(s) ~= sig_diffs_y(s)
            error('look at this')
        end
        dim = [start_x_pos+(x_width_per_mod*(s-1)) 0 .1 .1];
        str = x_labs2{s};
        annotation('textbox',dim,'String',str,'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    end
    
    legend([pym, pyf, pem, pef], 'young male', 'young female', 'elderly male', 'elderly female', 'Location', 'best')
    
    % save plot
    filename = ['fig4_ye_' overall_sig{1,1}];
    PrintFigs(gcf, paper_size(3:4)/100, [up.pub.save_folder, filename])
    
end

end

function compare_rr_mods(win_data, up)

fprintf('\n--- Fig 5');

% extract relevant index measures
eval(['rel_meas = win_data.' up.pub.rel_index ';']);

% cycle through each sig
for overall_sig = {'ppg', 'ekg'}
    
    % specify signal(s) of interest
    if strcmp(overall_sig, 'ppg')
        sig = 'ppgclin';
    elseif strcmp(overall_sig, 'ekg')
        sig = 'ekgclin';
    end
    
    % identify elements which are relevant to this signal
    eval(['rel_sig_els = win_data.' sig '_log & ~isnan(rel_meas);']);
    
    % identify rr cutoffs
    pops = {'rr1', 'rr2', 'rr3', 'rr4', 'rr5'};
    
    no_quantiles = length(pops);
    for pop_no = 1 : length(pops)
        quantile_no = str2double(pops{pop_no}(3));
        quantile_cutoffs.l = (quantile_no-1)/no_quantiles;
        quantile_cutoffs.u = quantile_no/no_quantiles;
        rr_cutoffs.l(pop_no) = quantile(win_data.ref(rel_sig_els), quantile_cutoffs.l);
        rr_cutoffs.u(pop_no) = quantile(win_data.ref(rel_sig_els), quantile_cutoffs.u);        
    end
    rr_cutoffs.u(length(pops)) = rr_cutoffs.u(length(pops))+0.000001;
    
    counter_no = 0; res = nan(1000,1000); [groups.mod, groups.age] = deal(cell(0)); 
    for pop_no = 1:length(pops)
        
        rel_pop_els = rel_sig_els & win_data.ref>= rr_cutoffs.l(pop_no) & win_data.ref< rr_cutoffs.u(pop_no);
        
        % specify modulations
        xb_mods = unique(win_data.m_xb(rel_sig_els)); xb_mods = xb_mods(~isnan(xb_mods));
        xa_mods = unique(win_data.m_xa(rel_sig_els)); xa_mods = xa_mods(~isnan(xa_mods));
        no_mods = length(xa_mods)+length(xb_mods);
        x_labs = cell(no_mods,1);
        for mod_no = 1 : no_mods
            if mod_no <= length(xa_mods)
                x_labs{mod_no} = ['XA' num2str(xa_mods(mod_no))];
            else
                x_labs{mod_no} = ['XB' num2str(xb_mods(mod_no-length(xa_mods)))];
            end
        end
        poss_mods = [xa_mods; xb_mods];
        
        for mod_no = 1 : length(poss_mods)
            curr_mod = poss_mods(mod_no);
            % identify relevant elements for this mod            
            if mod_no <= length(xa_mods)
                rel_mod_els = rel_sig_els & win_data.m_xa == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            else
                rel_mod_els = rel_sig_els & win_data.m_xb == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            end
            % identify index measures for these elements
            temp_ind_meas = rel_meas(rel_mod_els);
            % identify win nos for these elements
            temp_win_no = win_data.win_no(rel_mod_els);
            % look at each subject
            temp_subj = win_data.subj(rel_mod_els);
            subjs = unique(temp_subj);
            subj_med_ind_meas = nan(1,length(subjs));
            for subj_no = 1:length(subjs)
                curr_subj = subjs(subj_no);
                subj_wins = unique(temp_win_no(temp_subj == curr_subj));
                temp_meas = nan(length(subj_wins),1);
                for win_no = 1 : length(subj_wins)
                    curr_win = subj_wins(win_no);
                    cand_ind_meas = unique(temp_ind_meas(temp_win_no == curr_win & temp_subj == curr_subj));
                    if length(cand_ind_meas) == 1
                        temp_meas(win_no,1) = cand_ind_meas;
                    else
                        error('Check this')
                    end
                end
                subj_med_ind_meas(subj_no) = median(temp_meas);
            end
            counter_no = counter_no+1;
            groups.age{counter_no,1} = pops{pop_no};
            groups.mod{counter_no,1} = ['    $' x_labs{mod_no}(1), '_{', x_labs{mod_no}(2:end), '}$'];   % for the latex interpreter
            res(counter_no,1:length(subj_med_ind_meas)) = subj_med_ind_meas;
            clear subj_med_ind_meas subj_no curr_subj subjs temp_subj
            clear curr_mod rel_mod_els
        end
        clear mod_no temp_ind_meas
        
    end
    res = res(1:length(groups.age),:); res = res';
    
    % adjust res for young median values
    y_medians = nan(length(poss_mods),1);
    for mod_no = 1 : (length(poss_mods))
        young_mod_els = ~cellfun(@isempty, strfind(groups.age, 'rr3')) & strcmp(groups.mod, groups.mod{mod_no});
        temp_data = res(:,young_mod_els);
        y_medians(mod_no) = nanmedian(temp_data(:));        
    end
    for col_no = 1 : size(res, 2)
        mod_no = find(strcmp(groups.mod, groups.mod{col_no}));
        res(:,col_no) = res(:,col_no) - y_medians(mod_no(1));
    end
    
    % statistical analysis
    [mod_stats.p, mod_stats.z, mod_stats.tau] = deal(nan(length(poss_mods),1));
    for mod_no = 1 : (length(poss_mods))
        curr_mod = poss_mods(mod_no);
        % identify relevant elements for this mod
        if mod_no <= length(xa_mods)
            rel_mod_els = rel_sig_els & win_data.m_xa == curr_mod & ~isnan(rel_meas) & win_data.comb_log;
        else
            rel_mod_els = rel_sig_els & win_data.m_xb == curr_mod & ~isnan(rel_meas) & win_data.comb_log;
        end
        rel_alg_nos = unique(win_data.alg_no(rel_mod_els));
        rel_mod_els = rel_mod_els & win_data.alg_no == rel_alg_nos(1);
        rr_data = win_data.ref(rel_mod_els);
        meas_data = rel_meas(rel_mod_els);
        stats = mann_kendall_test(rr_data, meas_data);
        mod_stats.p(mod_no) = stats.p;
        mod_stats.z(mod_no) = stats.zval;
        mod_stats.tau(mod_no) = stats.tau;
    end
    % correction for multiple comparisons
    sig_diffs = correct_multiple_comparisons(mod_stats, up);
    
    % setup figure
    lwidth = 1;
    lbl_ftsize = 14;
    paper_size = [200, 200, 1000, 600];
    figure('Position',paper_size);
    
    % subplot 1
    subplot('Position', [0.22, 0.2, 0.75, 0.75])
    
    % make boxplot
    ylims = [-0.5, 0.5];
    new_labels = groups.mod(1:length(poss_mods));
    new_labels = [repmat({''}, [2,length(new_labels)]); new_labels'; repmat({''}, [2,length(new_labels)])]; new_labels = new_labels(:);
    boxplot(res,{groups.mod,groups.age},'colors','k','factorgap',[8 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
    hold on
    % add x-axis
    hline = refline([0 0]);
    hline.Color = 'k';
    % colour in boxes
    h = findobj(gca,'Tag','Box');
    box_diffs = repmat(sig_diffs(:)', [5,1]); sig_diffs = sig_diffs(:);
    box_ages = repmat(pops, [1, length(poss_mods)]);
    for j=1:length(h)
        if strcmp(box_ages(length(h)+1-j), 'rr1')
            prr1 = patch(get(h(j),'XData'),get(h(j),'YData'),0.9*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr2')
            prr2 = patch(get(h(j),'XData'),get(h(j),'YData'),0.75*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr3')
            prr3 = patch(get(h(j),'XData'),get(h(j),'YData'),0.6*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr4')
            prr4 = patch(get(h(j),'XData'),get(h(j),'YData'),0.45*[1,1,1]); % dark grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr5')
            prr5 = patch(get(h(j),'XData'),get(h(j),'YData'),0.3*[1,1,1]); % dark grey
        end
    end
    % replot boxplot
    boxplot(res,{groups.mod,groups.age},'colors','k','factorgap',[8 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
    for j=1:length(h)
        if box_diffs(length(h)+1-j) == 1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',[0.5,0.5,1], 'LineWidth',lwidth)
        elseif box_diffs(length(h)+1-j) == -1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',0.8*[1,1,0], 'LineWidth',lwidth)
        end
    end
    % colour in median line
    h = findobj(gca,'tag','Median');
    set(h,'Color','r')
    % y-axis label
    ylab = ylabel({'CC', 'relative to', 'mid quintile', 'median CC'}, 'FontSize', lbl_ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0], 'VerticalAlignment', 'middle');
    set(gca, 'FontSize', lbl_ftsize)
    box off
    grid on
    set(gca,'XGrid', 'off', 'TickLabelInterpreter', 'latex')
    ylim(up.pub.plot_lims), set(gca,'YTick', up.pub.y_ticks)
    
    % subplot 2
    left_hand_x_pos = 0.195; x_width = 0.77;
    annotation('textbox',[0 0.025 left_hand_x_pos .1],'String',{'Kendall''s', 'rank CC:'},'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    x_width_per_mod = x_width/length(poss_mods);
    start_x_pos = left_hand_x_pos - 0.002;
    x_labs2 = cell(1,length(poss_mods));
    for s = 1 : length(poss_mods)
        if sig_diffs(s) == 0
            x_labs2{s} = '';    % NS
        elseif sig_diffs(s) == 1
            x_labs2{s} = ['+' sprintf('%.2f', mod_stats.tau(s))];
        elseif sig_diffs(s) == -1
            x_labs2{s} = sprintf('%.2f', mod_stats.tau(s));
        end
        dim = [start_x_pos+(x_width_per_mod*(s-1)) 0 .1 .1];
        str = x_labs2{s};
        annotation('textbox',dim,'String',str,'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    end
    
    % create legend
    leg_labels = cell(0);
    for pop_no = 1 : length(pops)
        temp_txt = [ sprintf('%.0f', rr_cutoffs.l(pop_no)) ' \leq RR < ' sprintf('%.0f', rr_cutoffs.u(pop_no))];
        if length(temp_txt) == 14
            temp_txt = ['  ' temp_txt];
        end
        leg_labels{pop_no,1} = temp_txt;
    end
    leg_h = legend([prr1, prr2, prr3, prr4, prr5], leg_labels{1}, leg_labels{2}, leg_labels{3}, leg_labels{4}, leg_labels{5});
    set(leg_h, 'FontSize', lbl_ftsize - 4, 'Position',[0.01,0.8,0.16,0.17])
    
    % save plot
    filename = ['fig5_rr_' overall_sig{1,1}];
    PrintFigs(gcf, paper_size(3:4)/100, [up.pub.save_folder, filename])
    
end

end

function compare_rr_mods_ind_subjs(win_data, up)

fprintf('\n--- Fig 5');

% extract relevant index measures
eval(['rel_meas = win_data.' up.pub.rel_index ';']);

% cycle through each sig
for overall_sig = {'ppg', 'ekg'}
    
    % specify signal(s) of interest
    if strcmp(overall_sig, 'ppg')
        sig = 'ppgclin';
    elseif strcmp(overall_sig, 'ekg')
        sig = 'ekgclin';
    end
    
    % identify elements which are relevant to this signal
    eval(['rel_sig_els = win_data.' sig '_log & ~isnan(rel_meas);']);
    
    % identify rr cutoffs
    pops = {'rr1', 'rr2', 'rr3', 'rr4', 'rr5'};
    
    no_quantiles = length(pops);
    for pop_no = 1 : length(pops)
        quantile_no = str2double(pops{pop_no}(3));
        quantile_cutoffs.l = (quantile_no-1)/no_quantiles;
        quantile_cutoffs.u = quantile_no/no_quantiles;
        rr_cutoffs.l(pop_no) = quantile(win_data.ref(rel_sig_els), quantile_cutoffs.l);
        rr_cutoffs.u(pop_no) = quantile(win_data.ref(rel_sig_els), quantile_cutoffs.u);        
    end
    rr_cutoffs.u(length(pops)) = rr_cutoffs.u(length(pops))+0.000001;
    
    counter_no = 0; res = nan(1000,1000); [groups.mod, groups.age] = deal(cell(0)); 
    for pop_no = 1:length(pops)
        
        rel_pop_els = rel_sig_els & win_data.ref>= rr_cutoffs.l(pop_no) & win_data.ref< rr_cutoffs.u(pop_no);
        
        % specify modulations
        xb_mods = unique(win_data.m_xb(rel_sig_els)); xb_mods = xb_mods(~isnan(xb_mods));
        xa_mods = unique(win_data.m_xa(rel_sig_els)); xa_mods = xa_mods(~isnan(xa_mods));
        no_mods = length(xa_mods)+length(xb_mods);
        x_labs = cell(no_mods,1);
        for mod_no = 1 : no_mods
            if mod_no <= length(xa_mods)
                x_labs{mod_no} = ['Xa' num2str(xa_mods(mod_no))];
            else
                x_labs{mod_no} = ['Xb' num2str(xb_mods(mod_no-length(xa_mods)))];
            end
        end
        poss_mods = [xa_mods; xb_mods];
        
        for mod_no = 1 : length(poss_mods)
            curr_mod = poss_mods(mod_no);
            % identify relevant elements for this mod            
            if mod_no <= length(xa_mods)
                rel_mod_els = rel_sig_els & win_data.m_xa == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            else
                rel_mod_els = rel_sig_els & win_data.m_xb == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            end
            % identify index measures for these elements
            temp_ind_meas = rel_meas(rel_mod_els);
            % identify win nos for these elements
            temp_win_no = win_data.win_no(rel_mod_els);
            % look at each subject
            temp_subj = win_data.subj(rel_mod_els);
            subjs = unique(temp_subj);
            subj_med_ind_meas = nan(1,length(subjs));
            for subj_no = 1:length(subjs)
                curr_subj = subjs(subj_no);
                subj_wins = unique(temp_win_no(temp_subj == curr_subj));
                temp_meas = nan(length(subj_wins),1);
                for win_no = 1 : length(subj_wins)
                    curr_win = subj_wins(win_no);
                    cand_ind_meas = unique(temp_ind_meas(temp_win_no == curr_win & temp_subj == curr_subj));
                    if length(cand_ind_meas) == 1
                        temp_meas(win_no,1) = cand_ind_meas;
                    else
                        error('Check this')
                    end
                end
                subj_med_ind_meas(subj_no) = median(temp_meas);
            end
            counter_no = counter_no+1;
            groups.age{counter_no,1} = pops{pop_no};
            groups.mod{counter_no,1} = ['    ' x_labs{mod_no}];
            res(counter_no,1:length(subj_med_ind_meas)) = subj_med_ind_meas;
            clear subj_med_ind_meas subj_no curr_subj subjs temp_subj
            clear curr_mod rel_mod_els
            
        end
        clear mod_no temp_ind_meas
        
    end
    res = res(1:length(groups.age),:); res = res';
    
    % adjust res for young median values
    y_medians = nan(length(poss_mods),1);
    for mod_no = 1 : (length(poss_mods))
        young_mod_els = ~cellfun(@isempty, strfind(groups.age, 'rr3')) & strcmp(groups.mod, groups.mod{mod_no});
        temp_data = res(:,young_mod_els);
        y_medians(mod_no) = nanmedian(temp_data(:));        
    end
    for col_no = 1 : size(res, 2)
        mod_no = find(strcmp(groups.mod, groups.mod{col_no}));
        res(:,col_no) = res(:,col_no) - y_medians(mod_no(1));
    end
    
    % statistical analysis
    [mod_stats.p, mod_stats.z, mod_stats.tau] = deal(nan(length(poss_mods),1));
    for mod_no = 1 : (length(poss_mods))
        curr_mod = poss_mods(mod_no);
        % identify relevant elements for this mod
        if mod_no <= length(xa_mods)
            rel_mod_els = rel_sig_els & win_data.m_xa == curr_mod & ~isnan(rel_meas) & win_data.comb_log;
        else
            rel_mod_els = rel_sig_els & win_data.m_xb == curr_mod & ~isnan(rel_meas) & win_data.comb_log;
        end
        rel_alg_nos = unique(win_data.alg_no(rel_mod_els));
        rel_mod_els = rel_mod_els & win_data.alg_no == rel_alg_nos(1);
        rr_data = win_data.ref(rel_mod_els);
        meas_data = rel_meas(rel_mod_els);
        stats = mann_kendall_test(rr_data, meas_data);
        mod_stats.p(mod_no) = stats.p;
        mod_stats.z(mod_no) = stats.zval;
        mod_stats.tau(mod_no) = stats.tau;
    end
    % correction for multiple comparisons
    sig_diffs = correct_multiple_comparisons(mod_stats, up);
    
    % Perform analysis for individual subjects
    subjs = unique(win_data.subj);
    for subj_no = 1 : length(subjs)
        curr_subj = subjs(subj_no);
        for mod_no = 1 : (length(poss_mods))
            curr_mod = poss_mods(mod_no);
            % identify relevant elements for this mod
            if mod_no <= length(xa_mods)
                rel_mod_els = win_data.subj == curr_subj & rel_sig_els & win_data.m_xa == curr_mod & ~isnan(rel_meas) & win_data.comb_log;
            else
                rel_mod_els = win_data.subj == curr_subj & rel_sig_els & win_data.m_xb == curr_mod & ~isnan(rel_meas) & win_data.comb_log;
            end
            if sum(rel_mod_els) == 0 
                continue
            end
            rel_alg_nos = unique(win_data.alg_no(rel_mod_els));
            rel_mod_els = rel_mod_els & win_data.alg_no == rel_alg_nos(1);
            rr_data = win_data.ref(rel_mod_els);
            meas_data = rel_meas(rel_mod_els);
            if length(rr_data) <= 1
                continue
            end
            stats = mann_kendall_test(rr_data, meas_data);
            all_mod_stats.n(mod_no, subj_no) = length(rr_data);
            all_mod_stats.r(mod_no, subj_no) = range(rr_data);
            all_mod_stats.p(mod_no,subj_no) = stats.p;
            all_mod_stats.z(mod_no,subj_no) = stats.zval;
            all_mod_stats.tau(mod_no,subj_no) = stats.tau;
        end
    end
    
    % to test for individual subjects:
    rel_rows = sig_diffs==-1;
    prop_reduced_at_higher_RRs = sum(sum(all_mod_stats.tau(rel_rows,:)<0 & all_mod_stats.p(rel_rows,:)<0.05))/numel(all_mod_stats.p(rel_rows,:));
    prop_increased_at_higher_RRs = sum(sum(all_mod_stats.tau(rel_rows,:)>0 & all_mod_stats.p(rel_rows,:)<0.05))/numel(all_mod_stats.p(rel_rows,:));
    
    clear all_mod_stats
    
    % setup figure
    lwidth = 1;
    lbl_ftsize = 14;
    paper_size = [200, 200, 1000, 600];
    figure('Position',paper_size);
    
    % subplot 1
    subplot('Position', [0.22, 0.2, 0.75, 0.75])
    
    % make boxplot
    ylims = [-0.5, 0.5];
    new_labels = groups.mod(1:length(poss_mods));
    new_labels = [new_labels'; repmat({''}, [4,length(new_labels)])]; new_labels = new_labels(:);
    boxplot(res,{groups.mod,groups.age},'colors','k','factorgap',[8 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
    hold on
    % add x-axis
    hline = refline([0 0]);
    hline.Color = 'k';
    % colour in boxes
    h = findobj(gca,'Tag','Box');
    box_diffs = repmat(sig_diffs(:)', [5,1]); sig_diffs = sig_diffs(:);
    box_ages = repmat(pops, [1, length(poss_mods)]);
    for j=1:length(h)
        if strcmp(box_ages(length(h)+1-j), 'rr1')
            prr1 = patch(get(h(j),'XData'),get(h(j),'YData'),0.9*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr2')
            prr2 = patch(get(h(j),'XData'),get(h(j),'YData'),0.75*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr3')
            prr3 = patch(get(h(j),'XData'),get(h(j),'YData'),0.6*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr4')
            prr4 = patch(get(h(j),'XData'),get(h(j),'YData'),0.45*[1,1,1]); % dark grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr5')
            prr5 = patch(get(h(j),'XData'),get(h(j),'YData'),0.3*[1,1,1]); % dark grey
        end
    end
    % replot boxplot
    boxplot(res,{groups.mod,groups.age},'colors','k','factorgap',[8 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
    for j=1:length(h)
        if box_diffs(length(h)+1-j) == 1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',[0.5,0.5,1], 'LineWidth',lwidth)
        elseif box_diffs(length(h)+1-j) == -1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',0.8*[1,1,0], 'LineWidth',lwidth)
        end
    end
    % colour in median line
    h = findobj(gca,'tag','Median');
    set(h,'Color','r')
    % y-axis label
    ylab = ylabel({'CC', 'relative to', 'mid quintile', 'median CC'}, 'FontSize', lbl_ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0], 'VerticalAlignment', 'middle');
    set(gca, 'FontSize', lbl_ftsize)
    box off
    grid on
    set(gca,'XGrid', 'off')
    ylim(up.pub.plot_lims), set(gca,'YTick', up.pub.y_ticks)
    
    % subplot 2
    left_hand_x_pos = 0.195; x_width = 0.77;
    annotation('textbox',[0 0.025 left_hand_x_pos .1],'String',{'Kendall''s', 'rank CC:'},'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    x_width_per_mod = x_width/length(poss_mods);
    start_x_pos = left_hand_x_pos - 0.002;
    x_labs2 = cell(1,length(poss_mods));
    for s = 1 : length(poss_mods)
        if sig_diffs(s) == 0
            x_labs2{s} = '';    % NS
        elseif sig_diffs(s) == 1
            x_labs2{s} = ['+' sprintf('%.2f', mod_stats.tau(s))];
        elseif sig_diffs(s) == -1
            x_labs2{s} = sprintf('%.2f', mod_stats.tau(s));
        end
        dim = [start_x_pos+(x_width_per_mod*(s-1)) 0 .1 .1];
        str = x_labs2{s};
        annotation('textbox',dim,'String',str,'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    end
    
    % create legend
    leg_labels = cell(0);
    for pop_no = 1 : length(pops)
        temp_txt = [ sprintf('%.0f', rr_cutoffs.l(pop_no)) ' \leq RR < ' sprintf('%.0f', rr_cutoffs.u(pop_no))];
        if length(temp_txt) == 14
            temp_txt = ['  ' temp_txt];
        end
        leg_labels{pop_no,1} = temp_txt;
    end
    leg_h = legend([prr1, prr2, prr3, prr4, prr5], leg_labels{1}, leg_labels{2}, leg_labels{3}, leg_labels{4}, leg_labels{5});
    set(leg_h, 'FontSize', lbl_ftsize - 4, 'Position',[0.01,0.8,0.16,0.17])
    
    % save plot
    filename = ['fig5_rr_' overall_sig{1,1}];
    PrintFigs(gcf, paper_size(3:4)/100, [up.pub.save_folder, filename])
    
end

end

function compare_hr_mods(win_data, up)

fprintf('\n--- Fig 5');

% extract relevant index measures
eval(['rel_meas = win_data.' up.pub.rel_index ';']);

% cycle through each sig
for overall_sig = {'ppg', 'ekg'}
    
    % specify signal(s) of interest
    if strcmp(overall_sig, 'ppg')
        sig = 'ppgclin';
    elseif strcmp(overall_sig, 'ekg')
        sig = 'ekgclin';
    end
    
    % identify elements which are relevant to this signal
    eval(['rel_sig_els = win_data.' sig '_log & ~isnan(rel_meas);']);
    
    % identify rr cutoffs
    pops = {'rr1', 'rr2', 'rr3', 'rr4', 'rr5'};
    
    no_quantiles = length(pops);
    for pop_no = 1 : length(pops)
        quantile_no = str2double(pops{pop_no}(3));
        quantile_cutoffs.l = (quantile_no-1)/no_quantiles;
        quantile_cutoffs.u = quantile_no/no_quantiles;
        rr_cutoffs.l(pop_no) = quantile(win_data.hr(rel_sig_els), quantile_cutoffs.l);
        rr_cutoffs.u(pop_no) = quantile(win_data.hr(rel_sig_els), quantile_cutoffs.u);        
    end
    rr_cutoffs.u(length(pops)) = rr_cutoffs.u(length(pops))+0.000001;
    
    counter_no = 0; res = nan(1000,1000); [groups.mod, groups.age] = deal(cell(0)); 
    for pop_no = 1:length(pops)
        
        rel_pop_els = rel_sig_els & win_data.hr>= rr_cutoffs.l(pop_no) & win_data.hr< rr_cutoffs.u(pop_no);
        
        % specify modulations
        xb_mods = unique(win_data.m_xb(rel_sig_els)); xb_mods = xb_mods(~isnan(xb_mods));
        xa_mods = unique(win_data.m_xa(rel_sig_els)); xa_mods = xa_mods(~isnan(xa_mods));
        no_mods = length(xa_mods)+length(xb_mods);
        x_labs = cell(no_mods,1);
        for mod_no = 1 : no_mods
            if mod_no <= length(xa_mods)
                x_labs{mod_no} = ['Xa' num2str(xa_mods(mod_no))];
            else
                x_labs{mod_no} = ['Xb' num2str(xb_mods(mod_no-length(xa_mods)))];
            end
        end
        poss_mods = [xa_mods; xb_mods];
        
        for mod_no = 1 : length(poss_mods)
            curr_mod = poss_mods(mod_no);
            % identify relevant elements for this mod            
            if mod_no <= length(xa_mods)
                rel_mod_els = rel_sig_els & win_data.m_xa == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            else
                rel_mod_els = rel_sig_els & win_data.m_xb == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            end
            % identify index measures for these elements
            temp_ind_meas = rel_meas(rel_mod_els);
            % identify win nos for these elements
            temp_win_no = win_data.win_no(rel_mod_els);
            % look at each subject
            temp_subj = win_data.subj(rel_mod_els);
            subjs = unique(temp_subj);
            subj_med_ind_meas = nan(1,length(subjs));
            for subj_no = 1:length(subjs)
                curr_subj = subjs(subj_no);
                subj_wins = unique(temp_win_no(temp_subj == curr_subj));
                temp_meas = nan(length(subj_wins),1);
                for win_no = 1 : length(subj_wins)
                    curr_win = subj_wins(win_no);
                    cand_ind_meas = unique(temp_ind_meas(temp_win_no == curr_win & temp_subj == curr_subj));
                    if length(cand_ind_meas) == 1
                        temp_meas(win_no,1) = cand_ind_meas;
                    else
                        error('Check this')
                    end
                end
                subj_med_ind_meas(subj_no) = median(temp_meas);
            end
            counter_no = counter_no+1;
            groups.age{counter_no,1} = pops{pop_no};
            groups.mod{counter_no,1} = ['    ' x_labs{mod_no}];
            res(counter_no,1:length(subj_med_ind_meas)) = subj_med_ind_meas;
            clear subj_med_ind_meas subj_no curr_subj subjs temp_subj
            clear curr_mod rel_mod_els
        end
        clear mod_no temp_ind_meas
        
    end
    res = res(1:length(groups.age),:); res = res';
    
    % adjust res for young median values
    y_medians = nan(length(poss_mods),1);
    for mod_no = 1 : (length(poss_mods))
        young_mod_els = ~cellfun(@isempty, strfind(groups.age, 'rr3')) & strcmp(groups.mod, groups.mod{mod_no});
        temp_data = res(:,young_mod_els);
        y_medians(mod_no) = nanmedian(temp_data(:));        
    end
    for col_no = 1 : size(res, 2)
        mod_no = find(strcmp(groups.mod, groups.mod{col_no}));
        res(:,col_no) = res(:,col_no) - y_medians(mod_no(1));
    end
    
    % statistical analysis
    [mod_stats.p, mod_stats.z, mod_stats.tau] = deal(nan(length(poss_mods),1));
    for mod_no = 1 : (length(poss_mods))
        curr_mod = poss_mods(mod_no);
        % identify relevant elements for this mod
        if mod_no <= length(xa_mods)
            rel_mod_els = rel_sig_els & win_data.m_xa == curr_mod & ~isnan(rel_meas) & win_data.comb_log;
        else
            rel_mod_els = rel_sig_els & win_data.m_xb == curr_mod & ~isnan(rel_meas) & win_data.comb_log;
        end
        rel_alg_nos = unique(win_data.alg_no(rel_mod_els));
        rel_mod_els = rel_mod_els & win_data.alg_no == rel_alg_nos(1);
        rr_data = win_data.hr(rel_mod_els);
        meas_data = rel_meas(rel_mod_els);
        stats = mann_kendall_test(rr_data, meas_data);
        mod_stats.p(mod_no) = stats.p;
        mod_stats.z(mod_no) = stats.zval;
        mod_stats.tau(mod_no) = stats.tau;
    end
    % correction for multiple comparisons
    sig_diffs = correct_multiple_comparisons(mod_stats, up);
    
    % setup figure
    lwidth = 1;
    lbl_ftsize = 14;
    paper_size = [200, 200, 1000, 600];
    figure('Position',paper_size);
    
    % subplot 1
    subplot('Position', [0.22, 0.2, 0.75, 0.75])
    
    % make boxplot
    ylims = [-0.5, 0.5];
    new_labels = groups.mod(1:length(poss_mods));
    new_labels = [new_labels'; repmat({''}, [4,length(new_labels)])]; new_labels = new_labels(:);
    boxplot(res,{groups.mod,groups.age},'colors','k','factorgap',[8 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
    hold on
    % add x-axis
    hline = refline([0 0]);
    hline.Color = 'k';
    % colour in boxes
    h = findobj(gca,'Tag','Box');
    box_diffs = repmat(sig_diffs(:)', [5,1]); sig_diffs = sig_diffs(:);
    box_ages = repmat(pops, [1, length(poss_mods)]);
    for j=1:length(h)
        if strcmp(box_ages(length(h)+1-j), 'rr1')
            prr1 = patch(get(h(j),'XData'),get(h(j),'YData'),0.9*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr2')
            prr2 = patch(get(h(j),'XData'),get(h(j),'YData'),0.75*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr3')
            prr3 = patch(get(h(j),'XData'),get(h(j),'YData'),0.6*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr4')
            prr4 = patch(get(h(j),'XData'),get(h(j),'YData'),0.45*[1,1,1]); % dark grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr5')
            prr5 = patch(get(h(j),'XData'),get(h(j),'YData'),0.3*[1,1,1]); % dark grey
        end
    end
    % replot boxplot
    boxplot(res,{groups.mod,groups.age},'colors','k','factorgap',[8 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
    for j=1:length(h)
        if box_diffs(length(h)+1-j) == 1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',[0.5,0.5,1], 'LineWidth',lwidth)
        elseif box_diffs(length(h)+1-j) == -1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',0.8*[1,1,0], 'LineWidth',lwidth)
        end
    end
    % colour in median line
    h = findobj(gca,'tag','Median');
    set(h,'Color','r')
    % y-axis label
    ylab = ylabel({'CC', 'relative to', 'mid quintile', 'median CC'}, 'FontSize', lbl_ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0], 'VerticalAlignment', 'middle');
    set(gca, 'FontSize', lbl_ftsize)
    box off
    grid on
    set(gca,'XGrid', 'off')
    ylim(up.pub.plot_lims), set(gca,'YTick', up.pub.y_ticks)
    
    % subplot 2
    left_hand_x_pos = 0.195; x_width = 0.77;
    annotation('textbox',[0 0.025 left_hand_x_pos .1],'String',{'Kendall''s', 'rank CC:'},'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    x_width_per_mod = x_width/length(poss_mods);
    start_x_pos = left_hand_x_pos - 0.002;
    x_labs2 = cell(1,length(poss_mods));
    for s = 1 : length(poss_mods)
        if sig_diffs(s) == 0
            x_labs2{s} = '';    % NS
        elseif sig_diffs(s) == 1
            x_labs2{s} = ['+' sprintf('%.2f', mod_stats.tau(s))];
        elseif sig_diffs(s) == -1
            x_labs2{s} = sprintf('%.2f', mod_stats.tau(s));
        end
        dim = [start_x_pos+(x_width_per_mod*(s-1)) 0 .1 .1];
        str = x_labs2{s};
        annotation('textbox',dim,'String',str,'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    end
    
    % create legend
    leg_labels = cell(0);
    for pop_no = 1 : length(pops)
        temp_txt = [ sprintf('%.0f', rr_cutoffs.l(pop_no)) ' \leq HR < ' sprintf('%.0f', rr_cutoffs.u(pop_no))];
        if length(temp_txt) == 14
            temp_txt = ['  ' temp_txt];
        end
        leg_labels{pop_no,1} = temp_txt;
    end
    leg_h = legend([prr1, prr2, prr3, prr4, prr5], leg_labels{1}, leg_labels{2}, leg_labels{3}, leg_labels{4}, leg_labels{5});
    set(leg_h, 'FontSize', lbl_ftsize - 4, 'Position',[0.01,0.8,0.16,0.17])
    
    % save plot
    filename = ['fig5_hr_' overall_sig{1,1}];
    PrintFigs(gcf, paper_size(3:4)/100, [up.pub.save_folder, filename])
    
end

end

function compare_hrrr_mods(win_data, up)

fprintf('\n--- Fig 5');

% extract relevant index measures
eval(['rel_meas = win_data.' up.pub.rel_index ';']);

% cycle through each sig
for overall_sig = {'ppg', 'ekg'}
    
    % specify signal(s) of interest
    if strcmp(overall_sig, 'ppg')
        sig = 'ppgclin';
    elseif strcmp(overall_sig, 'ekg')
        sig = 'ekgclin';
    end
    
    % identify elements which are relevant to this signal
    eval(['rel_sig_els = win_data.' sig '_log & ~isnan(rel_meas);']);
    
    % identify rr cutoffs
    pops = {'rr1', 'rr2', 'rr3', 'rr4', 'rr5'};
    
    no_quantiles = length(pops);
    for pop_no = 1 : length(pops)
        quantile_no = str2double(pops{pop_no}(3));
        quantile_cutoffs.l = (quantile_no-1)/no_quantiles;
        quantile_cutoffs.u = quantile_no/no_quantiles;
        rr_cutoffs.l(pop_no) = quantile(win_data.hr_rr(rel_sig_els), quantile_cutoffs.l);
        rr_cutoffs.u(pop_no) = quantile(win_data.hr_rr(rel_sig_els), quantile_cutoffs.u);        
    end
    rr_cutoffs.u(length(pops)) = rr_cutoffs.u(length(pops))+0.000001;
    
    counter_no = 0; res = nan(1000,1000); [groups.mod, groups.age] = deal(cell(0)); 
    for pop_no = 1:length(pops)
        
        rel_pop_els = rel_sig_els & win_data.hr_rr>= rr_cutoffs.l(pop_no) & win_data.hr_rr< rr_cutoffs.u(pop_no);
        
        % specify modulations
        xb_mods = unique(win_data.m_xb(rel_sig_els)); xb_mods = xb_mods(~isnan(xb_mods));
        xa_mods = unique(win_data.m_xa(rel_sig_els)); xa_mods = xa_mods(~isnan(xa_mods));
        no_mods = length(xa_mods)+length(xb_mods);
        x_labs = cell(no_mods,1);
        for mod_no = 1 : no_mods
            if mod_no <= length(xa_mods)
                x_labs{mod_no} = ['Xa' num2str(xa_mods(mod_no))];
            else
                x_labs{mod_no} = ['Xb' num2str(xb_mods(mod_no-length(xa_mods)))];
            end
        end
        poss_mods = [xa_mods; xb_mods];
        
        for mod_no = 1 : length(poss_mods)
            curr_mod = poss_mods(mod_no);
            % identify relevant elements for this mod            
            if mod_no <= length(xa_mods)
                rel_mod_els = rel_sig_els & win_data.m_xa == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            else
                rel_mod_els = rel_sig_els & win_data.m_xb == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            end
            % identify index measures for these elements
            temp_ind_meas = rel_meas(rel_mod_els);
            % identify win nos for these elements
            temp_win_no = win_data.win_no(rel_mod_els);
            % look at each subject
            temp_subj = win_data.subj(rel_mod_els);
            subjs = unique(temp_subj);
            subj_med_ind_meas = nan(1,length(subjs));
            for subj_no = 1:length(subjs)
                curr_subj = subjs(subj_no);
                subj_wins = unique(temp_win_no(temp_subj == curr_subj));
                temp_meas = nan(length(subj_wins),1);
                for win_no = 1 : length(subj_wins)
                    curr_win = subj_wins(win_no);
                    cand_ind_meas = unique(temp_ind_meas(temp_win_no == curr_win & temp_subj == curr_subj));
                    if length(cand_ind_meas) == 1
                        temp_meas(win_no,1) = cand_ind_meas;
                    else
                        error('Check this')
                    end
                end
                subj_med_ind_meas(subj_no) = median(temp_meas);
            end
            counter_no = counter_no+1;
            groups.age{counter_no,1} = pops{pop_no};
            groups.mod{counter_no,1} = ['    ' x_labs{mod_no}];
            res(counter_no,1:length(subj_med_ind_meas)) = subj_med_ind_meas;
            clear subj_med_ind_meas subj_no curr_subj subjs temp_subj
            clear curr_mod rel_mod_els
        end
        clear mod_no temp_ind_meas
        
    end
    res = res(1:length(groups.age),:); res = res';
    
    % adjust res for young median values
    y_medians = nan(length(poss_mods),1);
    for mod_no = 1 : (length(poss_mods))
        young_mod_els = ~cellfun(@isempty, strfind(groups.age, 'rr3')) & strcmp(groups.mod, groups.mod{mod_no});
        temp_data = res(:,young_mod_els);
        y_medians(mod_no) = nanmedian(temp_data(:));        
    end
    for col_no = 1 : size(res, 2)
        mod_no = find(strcmp(groups.mod, groups.mod{col_no}));
        res(:,col_no) = res(:,col_no) - y_medians(mod_no(1));
    end
    
    % statistical analysis
    [mod_stats.p, mod_stats.z, mod_stats.tau] = deal(nan(length(poss_mods),1));
    for mod_no = 1 : (length(poss_mods))
        curr_mod = poss_mods(mod_no);
        % identify relevant elements for this mod
        if mod_no <= length(xa_mods)
            rel_mod_els = rel_sig_els & win_data.m_xa == curr_mod & ~isnan(rel_meas) & win_data.comb_log;
        else
            rel_mod_els = rel_sig_els & win_data.m_xb == curr_mod & ~isnan(rel_meas) & win_data.comb_log;
        end
        rel_alg_nos = unique(win_data.alg_no(rel_mod_els));
        rel_mod_els = rel_mod_els & win_data.alg_no == rel_alg_nos(1);
        rr_data = win_data.hr_rr(rel_mod_els);
        meas_data = rel_meas(rel_mod_els);
        stats = mann_kendall_test(rr_data, meas_data);
        mod_stats.p(mod_no) = stats.p;
        mod_stats.z(mod_no) = stats.zval;
        mod_stats.tau(mod_no) = stats.tau;
    end
    % correction for multiple comparisons
    sig_diffs = correct_multiple_comparisons(mod_stats, up);
    
    % setup figure
    lwidth = 1;
    lbl_ftsize = 14;
    paper_size = [200, 200, 1000, 600];
    figure('Position',paper_size);
    
    % subplot 1
    subplot('Position', [0.22, 0.2, 0.75, 0.75])
    
    % make boxplot
    ylims = [-0.5, 0.5];
    new_labels = groups.mod(1:length(poss_mods));
    new_labels = [new_labels'; repmat({''}, [4,length(new_labels)])]; new_labels = new_labels(:);
    boxplot(res,{groups.mod,groups.age},'colors','k','factorgap',[8 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
    hold on
    % add x-axis
    hline = refline([0 0]);
    hline.Color = 'k';
    % colour in boxes
    h = findobj(gca,'Tag','Box');
    box_diffs = repmat(sig_diffs(:)', [5,1]); sig_diffs = sig_diffs(:);
    box_ages = repmat(pops, [1, length(poss_mods)]);
    for j=1:length(h)
        if strcmp(box_ages(length(h)+1-j), 'rr1')
            prr1 = patch(get(h(j),'XData'),get(h(j),'YData'),0.9*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr2')
            prr2 = patch(get(h(j),'XData'),get(h(j),'YData'),0.75*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr3')
            prr3 = patch(get(h(j),'XData'),get(h(j),'YData'),0.6*[1,1,1]); % light grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr4')
            prr4 = patch(get(h(j),'XData'),get(h(j),'YData'),0.45*[1,1,1]); % dark grey
        elseif strcmp(box_ages(length(h)+1-j), 'rr5')
            prr5 = patch(get(h(j),'XData'),get(h(j),'YData'),0.3*[1,1,1]); % dark grey
        end
    end
    % replot boxplot
    boxplot(res,{groups.mod,groups.age},'colors','k','factorgap',[8 0],'labelverbosity','major', 'Widths', 0.8, 'Labels', new_labels); ylim(ylims)
    for j=1:length(h)
        if box_diffs(length(h)+1-j) == 1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',[0.5,0.5,1], 'LineWidth',lwidth)
        elseif box_diffs(length(h)+1-j) == -1
            lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
            xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
            ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
            rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',0.8*[1,1,0], 'LineWidth',lwidth)
        end
    end
    % colour in median line
    h = findobj(gca,'tag','Median');
    set(h,'Color','r')
    % y-axis label
    ylab = ylabel({'CC', 'relative to', 'mid quintile', 'median CC'}, 'FontSize', lbl_ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0], 'VerticalAlignment', 'middle');
    set(gca, 'FontSize', lbl_ftsize)
    box off
    grid on
    set(gca,'XGrid', 'off')
    ylim(up.pub.plot_lims), set(gca,'YTick', up.pub.y_ticks)
    
    % subplot 2
    left_hand_x_pos = 0.195; x_width = 0.77;
    annotation('textbox',[0 0.025 left_hand_x_pos .1],'String',{'Kendall''s', 'rank CC:'},'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    x_width_per_mod = x_width/length(poss_mods);
    start_x_pos = left_hand_x_pos - 0.002;
    x_labs2 = cell(1,length(poss_mods));
    for s = 1 : length(poss_mods)
        if sig_diffs(s) == 0
            x_labs2{s} = '';    % NS
        elseif sig_diffs(s) == 1
            x_labs2{s} = ['+' sprintf('%.2f', mod_stats.tau(s))];
        elseif sig_diffs(s) == -1
            x_labs2{s} = sprintf('%.2f', mod_stats.tau(s));
        end
        dim = [start_x_pos+(x_width_per_mod*(s-1)) 0 .1 .1];
        str = x_labs2{s};
        annotation('textbox',dim,'String',str,'LineStyle','none', 'FontSize', lbl_ftsize, 'HorizontalAlignment', 'center')
    end
    
    % create legend
    leg_labels = cell(0);
    for pop_no = 1 : length(pops)
        temp_txt = [ sprintf('%.0f', rr_cutoffs.l(pop_no)) ' \leq HR:RR < ' sprintf('%.0f', rr_cutoffs.u(pop_no))];
        if length(temp_txt) == 14
            temp_txt = ['  ' temp_txt];
        end
        leg_labels{pop_no,1} = temp_txt;
    end
    leg_h = legend([prr1, prr2, prr3, prr4, prr5], leg_labels{1}, leg_labels{2}, leg_labels{3}, leg_labels{4}, leg_labels{5});
    set(leg_h, 'FontSize', lbl_ftsize - 4, 'Position',[0.01,0.8,0.18,0.17])
    
    % save plot
    filename = ['fig5_hrrr_' overall_sig{1,1}];
    PrintFigs(gcf, paper_size(3:4)/100, [up.pub.save_folder, filename])
    
end

end

function compare_ecg_and_ppg_mods(win_data, up)

fprintf('\n--- Fig 6');

% specify signals of interest
if strcmp(up.pub.equip_type, 'clin')
    sigs = {'ppgclin', 'ekgclin'};
elseif strcmp(up.pub.equip_type, 'lab')
    sigs = {'ppgfraw', 'ekgraw'};
end

% extract relevant index measures
eval(['rel_meas = win_data.' up.pub.rel_index ';']);

pops = {'all', 'young', 'elderly'};
for pop_no = 1:length(pops)
    
    
    switch pops{pop_no}
        case 'all'
            rel_pop_els = true(length(win_data.subj),1);
        case 'young'
            rel_pop_els = win_data.young_log;
        case 'elderly'
            rel_pop_els = win_data.elderly_log;
    end
    
    % cycle through signals
    counter_no = 0; stats.med = nan(1000,1); [stats.mod, stats.sig] = deal(cell(0));
    for sig_no = 1 : length(sigs)
        
        % identify elements which are relevant to this signal
        eval(['rel_sig_els = win_data.' sigs{sig_no} '_log;']);
        
        % specify modulations
        xb_mods = unique(win_data.m_xb(rel_sig_els)); xb_mods = xb_mods(~isnan(xb_mods));
        xa_mods = unique(win_data.m_xa(rel_sig_els)); xa_mods = xa_mods(~isnan(xa_mods));
        no_mods = length(xa_mods)+length(xb_mods);
        x_labs = cell(no_mods,1);
        for mod_no = 1 : no_mods
            if mod_no <= length(xa_mods)
                x_labs{mod_no} = ['XA' num2str(xa_mods(mod_no))];
            else
                x_labs{mod_no} = ['XB' num2str(xb_mods(mod_no-length(xa_mods)))];
            end
        end
        poss_mods = [xa_mods; xb_mods];
        
        for mod_no = 1 : length(poss_mods)
            curr_mod = poss_mods(mod_no);
            % identify relevant elements for this mod            
            if mod_no <= length(xa_mods)
                rel_mod_els = rel_sig_els & win_data.m_xa == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            else
                rel_mod_els = rel_sig_els & win_data.m_xb == curr_mod & ~isnan(rel_meas) & win_data.comb_log & rel_pop_els;
            end
            % identify index measures for these elements
            temp_ind_meas = rel_meas(rel_mod_els);
            % identify win nos for these elements
            temp_win_no = win_data.win_no(rel_mod_els);
            % look at each subject
            temp_subj = win_data.subj(rel_mod_els);
            subjs = unique(temp_subj);
            subj_med_ind_meas = nan(1,max(subjs));
            for subj_no = 1:length(subjs)
                curr_subj = subjs(subj_no);
                subj_wins = unique(temp_win_no(temp_subj == curr_subj));
                temp_meas = nan(length(subj_wins),1);
                for win_no = 1 : length(subj_wins)
                    curr_win = subj_wins(win_no);
                    cand_ind_meas = unique(temp_ind_meas(temp_win_no == curr_win & temp_subj == curr_subj));
                    if length(cand_ind_meas) == 1
                        temp_meas(win_no,1) = cand_ind_meas;
                    else
                        error('Check this')
                    end
                end
                subj_med_ind_meas(subjs(subj_no)) = median(temp_meas);
            end
            counter_no = counter_no+1;
            stats.med(counter_no,1:length(subj_med_ind_meas)) = subj_med_ind_meas;
            stats.mod{counter_no} = ['    $' x_labs{mod_no}(1), '_{', x_labs{mod_no}(2:end), '}$'];   % for the latex interpreter
            stats.sig{counter_no} = sigs{sig_no}(1:3);
            clear subj_med_ind_meas subj_no curr_subj subjs temp_subj curr_mod rel_mod_els
        end
        clear mod_no rel_sig_els temp_ind_meas
        
    end
    clear sig_no
    stats.med = stats.med(1:counter_no,:);
    stats.med = stats.med';
    
    % eliminate any subjects with nans (which means they didn't contribute
    % a value for at least one mod)
    keep_log = true(size(stats.med,1),1);
    for subj_no = 1 : size(stats.med,1)
        rel_data = stats.med(subj_no,:);
        if sum(isnan(rel_data))
            keep_log(subj_no) = false;
        end
    end
    stats.med = stats.med(keep_log,:);
    
    % statistical analysis
    meds = median(stats.med); [~, rel_el] = max(meds);
    control_data = stats.med(:,rel_el); mod_stats.p = nan(size(stats.med,2),1);
    for mod_no = 1 : size(stats.med,2)
        comparison_data = stats.med(:,mod_no);
        test_data = [control_data(:); comparison_data(:)];
        test_groups = [repmat({'control'}, [length(control_data) 1]); repmat({'comparison'}, [length(comparison_data) 1])];
        mod_stats.p(mod_no) = kruskalwallis(test_data,test_groups,'off');
    end
    
    % order according to medians
    [~, order] = sort(mod_stats.p, 'descend');
    stats.med = stats.med(:,order);
    stats.mod = stats.mod(order);
    stats.sig = stats.sig(order);
    stats.sig = strrep(stats.sig, 'ekg' , 'ECG');
    stats.sig = strrep(stats.sig, 'ppg' , 'PPG');
    mod_stats.p = mod_stats.p(order);
    mod_stats.z = [1, -1*ones(1, size(stats.med,2)-1)];
    
    
    % correct for multiple comparisons
    sig_diffs = correct_multiple_comparisons(mod_stats, up);
    sig_diffs(1)=1;
    
    % setup figure
    lbl_ftsize = 14;
    paper_size = [200, 200, 1000, 500];
    figure('Position',paper_size);
    
    % subplot 1
    subplot('Position', [0.2, 0.2, 0.75, 0.75])
    
    % boxplot
    boxplot(stats.med), hold on
    % colour in boxes
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        if strcmp(stats.sig{length(h)+1-j}, 'PPG')
            pppg = patch(get(h(j),'XData'),get(h(j),'YData'),0.8*[1,1,1]); % light grey
        elseif strcmp(stats.sig{length(h)+1-j}, 'ECG')
            pecg = patch(get(h(j),'XData'),get(h(j),'YData'),0.4*[1,1,1]); % dark grey
        end
    end
    % y-axis label
    ylab = ylabel('CC', 'FontSize', lbl_ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
    set(gca, 'FontSize', lbl_ftsize)
    boxplot(stats.med, 'colors','k')
    
%     % outline statistical differences
%     lwidth = 2;
%     for j=1:length(h)
%         if sig_diffs(length(h)+1-j) == 1
%             lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
%             xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
%             ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
%             rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',[0.5,0.5,1], 'LineWidth',lwidth+1)
%         elseif sig_diffs(length(h)+1-j) == -1
%             lh_corner = [min(get(h(j),'XData')), min(get(h(j),'YData'))];
%             xlen = max(get(h(j),'XData')) - min(get(h(j),'XData'));
%             ylen = max(get(h(j),'YData')) - min(get(h(j),'YData'));
%             rectangle('Position', [lh_corner, xlen, ylen], 'EdgeColor',0.8*[1,1,0], 'LineWidth',lwidth)
%         end
%     end
    
    % create annotation of statistical differences
    sig_diff_els = find(sig_diffs == -1);
    first_sig_diff = sig_diff_els(1);
    last_sig_diff = sig_diff_els(end);
    if length(sig_diff_els) ~= (last_sig_diff-first_sig_diff+1)
        error('Check this')
    end
    start_x = min(get(h(length(h)-first_sig_diff+1),'XData'));
    end_x = max(get(h(length(h)-last_sig_diff+1),'XData'));
    temp=get(gca,'position');
    axlims=axis;
    start_coord_x = temp(1) + temp(3)*((start_x-axlims(1))/(axlims(2)-axlims(1)));
    end_coord_x = temp(1) + temp(3)*((end_x-axlims(1))/(axlims(2)-axlims(1)));
    start_y = 0.06;
    annotation('line', [start_coord_x end_coord_x], [start_y , start_y ])
    annotation('line', [start_coord_x start_coord_x], [start_y , start_y+0.04])
    annotation('line', [end_coord_x end_coord_x], [start_y , start_y+0.04])
    mid_coord_x = mean([start_coord_x, end_coord_x])-0.12;
    dim = [mid_coord_x 0.01 0.25 0.04];
    str = 'Significantly lower CCs';
    annotation('textbox',dim,'String',str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'LineStyle', 'none', 'FontSize', lbl_ftsize);
    
    % create annotation of control group
    dim = [0.085 0.01 0.25 0.04];
    str = 'Control';
    annotation('textbox',dim,'String',str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'LineStyle', 'none', 'FontSize', lbl_ftsize);
    annotation('line', [0.215 0.215], [start_y , start_y+0.04])
    
    % colour in median line
    h = findobj(gca,'tag','Median');
    set(h,'Color','r')
    
    % add mod labels
    set(gca,'XTickLabel',stats.mod, 'FontSize', lbl_ftsize, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'latex')
    box off
    grid on
    set(gca,'XGrid', 'off')
    ylim([0 1]), set(gca,'YTick', 0:0.2:1)
    
    % create legend
    legend([pecg, pppg], 'ECG', 'PPG', 'Location', 'best');
    
%     % Find out modulation types:
%     bw_log = ~cellfun(@isempty, strfind(stats.mod, 'Xa1')) ...
%         | ~cellfun(@isempty, strfind(stats.mod, 'Xb1'))...
%         | ~cellfun(@isempty, strfind(stats.mod, 'Xb4'))...
%         | ~cellfun(@isempty, strfind(stats.mod, 'Xb5'))...
%         | ~cellfun(@isempty, strfind(stats.mod, 'Xb6'))...
%         | ~cellfun(@isempty, strfind(stats.mod, 'Xb9'));
%     am_log = ~cellfun(@isempty, strfind(stats.mod, 'Xa2')) ...
%         | ~cellfun(@isempty, strfind(stats.mod, 'Xb2'))...
%         | ~cellfun(@isempty, strfind(stats.mod, 'Xb5'))...
%         | ~cellfun(@isempty, strfind(stats.mod, 'Xb6'))...
%         | ~cellfun(@isempty, strfind(stats.mod, 'Xb8'));
%     fm_log = ~cellfun(@isempty, strfind(stats.mod, 'Xa3')) ...
%         | ~cellfun(@isempty, strfind(stats.mod, 'Xb3'))...
%         | ~cellfun(@isempty, strfind(stats.mod, 'Xb7'))...
%         | ~cellfun(@isempty, strfind(stats.mod, 'Xb8'));
%     
%     % Annotate modulation types
%     start_x = min(mean(cell2mat(get(h,'XData')),2));
%     end_x = max(mean(cell2mat(get(h,'XData')),2));
%     temp=get(gca,'position');
%     x_dims = linspace(start_x, end_x, (length(h)/2));
%     axlims=axis;
%     
%     for x_dim_no = 1 : length(x_dims)
%         curr_x_dim = temp(1) - 0.12 + temp(3)*((x_dims(x_dim_no)-axlims(1))/(axlims(2)-axlims(1)));
%         % skip if more than one modulation
%         if bw_log(x_dim_no) + am_log(x_dim_no) + fm_log(x_dim_no) > 1
%             continue
%         end
%         % skip if this modulation has been done for this signal
%         if bw_log(x_dim_no)
%             rel_log = bw_log;
%         elseif am_log(x_dim_no)
%             rel_log = am_log;
%         elseif fm_log(x_dim_no)
%             rel_log = fm_log;
%         end
%         if sum(strcmp(stats.sig(1:(x_dim_no-1)), stats.sig{x_dim_no}) & ...
%                 rel_log(1:(x_dim_no-1)) == rel_log(x_dim_no))
%             continue
%         end
%         if bw_log(x_dim_no)
%             dim = [curr_x_dim 0.87 0.25 0.04];
%             str = 'BW';
%             annotation('textbox',dim,'String',str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'LineStyle', 'none', 'FontSize', lbl_ftsize);
%         end
%         if am_log(x_dim_no)
%             dim = [curr_x_dim 0.91 0.25 0.04];
%             str = 'AM';
%             annotation('textbox',dim,'String',str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'LineStyle', 'none', 'FontSize', lbl_ftsize);
%         end
%         if fm_log(x_dim_no)
%             dim = [curr_x_dim 0.95 0.25 0.04];
%             str = 'FM';
%             annotation('textbox',dim,'String',str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'LineStyle', 'none', 'FontSize', lbl_ftsize);
%         end
%         
%     end
    
    % save plot
    filename = ['fig6_' pops{pop_no}];
    PrintFigs(gcf, paper_size(3:4)/100, [up.pub.save_folder, filename])
    
    % clear vars
    clear x_vals med_vals lq_vals uq_vals x_labs ylab x_labs stats poss_mods
    
end

end

function stats = mann_kendall_test(pred, resp)

%% Mann-Kendall Test for Monotonic Trend
% designed using information at:
%    http://vsp.pnnl.gov/help/Vsample/Design_Trend_Mann_Kendall.htm

% sort data according to predictor variable
[pred, order] = sort(pred);
resp = resp(order);

% find s
s = 0;
for i = 1 : (length(resp)-1)
    for j = i: length(resp)
        s = s + sign(resp(j)-resp(i));
    end
end

% determine all possible differences between response variable measurements
resp1 = repmat(resp, [1,length(resp)]);
resp2 = repmat(resp', [length(resp),1]);
diffs = resp1-resp2;
rows = 1:size(diffs,1);
for col_no = 1 : size(diffs,2)
    exc_rows = rows<=col_no;
    diffs(exc_rows, col_no) = nan;
end
diffs = diffs(~isnan(diffs));

% allocate signs to differences
diff_signs = sign(diffs);

% find S: no of positive differences - no of negative differences
S = sum(diff_signs>0) - sum(diff_signs<0);

if S ~= s
    error('check this')
end

% assuming no of observations > 10, find variance of S
n = length(resp);   % no of observations
V_non_tied = (n*(n-1)*((2*n)+5)/18);
[resp_freqs,resp_vals]=hist(resp,unique(resp));
m = sum(resp_freqs>1);   % no of groups of tied ranks
t = resp_freqs(resp_freqs>1); % no of tied observations in each group of tied ranks
V_tied = sum(t.*(t-1).*((2*t)+5)/18);
V = V_non_tied + V_tied;

% find standardized variable u
if S > 0
    u = (S-1)/sqrt(V);
elseif S == 0
    u = 0;
elseif S < 0
    u = (S+1)/sqrt(V);
end

% find relevant threshold u-value for significance (using normal
% distribution and two-tailed test)
sig_level = 0.05;
P = 1-(sig_level/2);
u_thresh = norminv(P,0,1);
if abs(u) <= u_thresh
    h = 0;  % accept null hypothesis of no trend
elseif u> u_thresh
    h = 1;  % reject null hypothesis and accept alternative hypothesis of positive trend
elseif u < u_thresh
    h = -1; % reject null hypothesis and accept alternative hypothesis of negative trend
end

% find the Kendall rank correlation coefficient (tau):
max_score = length(resp)*(length(resp)-1)/2;
tau = s/max_score;

% output stats
stats.h = h;
stats.zval = u;
stats.p = 1-(abs(normcdf(abs(u))-normcdf(-1*abs(u))));
stats.tau = tau;

end

function sig_diffs = correct_multiple_comparisons(stats, up)

% extract required data, and order according to p-value (ascending)
[~,order] = sort(stats.p);
p = stats.p(order);
z = stats.z(order);

% setup variable
temp_sig_diffs = zeros(length(p),1);

% find significant differences
for comparison_no = 1 : length(p)
    % Holm's correction for only the number of remaining comparisons
    no_tests = length(p)-comparison_no+1;
    % Sidak's correction for the number of remaining comparisons
    alpha_sidak = 1 - ((1-up.pub.alpha)^(1/no_tests));
    if p(comparison_no) < alpha_sidak
        if z(comparison_no) < 0
            temp_sig_diffs(comparison_no) = -1;
        else
            temp_sig_diffs(comparison_no) = 1;
        end
    else
        break
    end
end

% return significant differences back to original order
orig_p(order) = p;
sig_diffs(order) = temp_sig_diffs;

end