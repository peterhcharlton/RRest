function specific_vortal_plots(up)

save_name = up.paths.filenames.hr_rr_scatter;
savepath = [up.paths.plots_save_folder, up.paths.filenames.hr_rr_scatter, '.eps'];
if ~up.analysis.redo_stats
    exist_log = check_exists(savepath, save_name);
    if exist_log
        return
    end
end

%% Make plot of feature vs filter resp sigs
if ~isempty(strfind(up.paths.paper_figures_folder, 'pc13')) & ...
        isempty(strfind(up.paths.data_save_folder, 'REC')) & isempty(strfind(up.paths.data_save_folder, 'WALK')) ...
        & isempty(strfind(up.paths.data_save_folder, 'EX'))
    % Load resp sigs
    subj = 2;
    load_path = [strrep(up.paths.data_save_folder, 'rest_and_rec', 'rest'), num2str(subj), up.paths.filenames.respSigs];
    resp_sigs = load(load_path);
    if sum(strcmp(fieldnames(resp_sigs), 'ppg_flt_Wam'))
        % Load raw sig
        load_path = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.int_respSigs];
        int_resp_sigs = load(load_path);
        % Load annotated breaths
        up.paths.root_ann_folder = 'C:\Documents\Data\VORTAL\Manual_breath_annotations\';
        up.paths.observers = {'TB', 'DV'};
        up.paths.ann_folder = '_annotations\';
        obs_no = 1;
        data_folder = [up.paths.root_ann_folder, up.paths.observers{obs_no}, up.paths.ann_folder];
        % find SID
        temp = num2str(subj);
        if length(temp) == 1
            SID = ['00', temp];
        elseif length(temp) == 2
            SID = ['0', temp];
        end
        % Load annotated breaths
        rel_name = ['VORTAL' SID 'rest' up.paths.observers{obs_no} '-paw_an'];
        anns = load([data_folder, rel_name]); t = anns.PKS.t;
        % Load s ref time
        up.paths.s_ref_time_folder = 'C:\Documents\Data\VORTAL\Analysis_files\Processed_Data\';
        loadpath = [up.paths.s_ref_time_folder, SID, 's_con_rest'];
        load(loadpath, 'Sreftime', 'Send', 'Sstart');
        load(up.paths.db_data, 'db_data');
        bst_time = db_data.bst_log(db_data.sid == subj)/(24);
        % find t, measured in seconds since the start of the period
        t = 60*60*24*(t/(60*60*24) - floor( t/(60*60*24) ) - Sreftime + bst_time) - Sstart; % As the times are in seconds since 1970
        % Identify relevant section
        start_t = 92;
        rel_t = [start_t, start_t+up.paramSet.winLeng];
        rel_feat = int_resp_sigs.ppg_FMeam_FPt_PDtIMS_EHF;
        rel_filt = resp_sigs.ppg_flt_Wam;
        rel_sig = int_resp_sigs.ppg_EHF;
        rel_els = rel_feat.t >= rel_t(1) & rel_feat.t <= rel_t(2);
        rel_data.feat.t = rel_feat.t(rel_els);
        rel_data.feat.v = rel_feat.v(rel_els);
        rel_els = rel_filt.t >= rel_t(1) & rel_filt.t <= rel_t(2);
        rel_data.filt.t = rel_filt.t(rel_els);
        rel_data.filt.v = rel_filt.v(rel_els);
        rel_els = rel_sig.t >= rel_t(1) & rel_sig.t <= rel_t(2);
        rel_data.sig.t = rel_sig.t(rel_els);
        rel_data.sig.v = rel_sig.v(rel_els);
        
        % Make Figure
        h_fig = figure('Position', [200, 200, 600, 300]);
        
        % Plot
        ftsize = 12; lwidth = 2;
        h1 = plot(rel_data.sig.t-rel_t(1), detrend(rel_data.sig.v), 'b', 'LineWidth',1);
        hold on
        plot(rel_data.filt.t-rel_t(1), detrend(rel_data.filt.v), 'k', 'LineWidth',lwidth+2)
        h2 = plot(rel_data.filt.t-rel_t(1), detrend(rel_data.filt.v), 'c', 'LineWidth',lwidth);
        h3 = plot(rel_data.feat.t-rel_t(1), detrend(rel_data.feat.v), '-r', 'LineWidth',lwidth);
        h4 = plot(t-start_t, -0.05+zeros(length(t),1), '.k','MarkerSize',30);
        plot(rel_data.feat.t-rel_t(1), detrend(rel_data.feat.v), '.r', 'LineWidth',lwidth,'MarkerSize',20)
        xlim(rel_t-rel_t(1))
        temp = range(rel_data.sig.v);
        ylim([min(rel_data.sig.v)-0.05*temp, max(rel_data.sig.v)+0.05*temp])
        xlabel('Time [s]', 'FontSize', ftsize)
        ylabel('PPG', 'FontSize', ftsize)
        legend([h1, h2, h3, h4], {'PPG', 'Filter', 'Feature', 'Breaths'},'Location','northoutside','Orientation','horizontal')
        set(gca, 'FontSize', ftsize)
        set(gca, 'YTick', [])
        % Save
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperSize', [6, 3]);
        set(gcf,'PaperPosition',[0 0 6 3]);
        save_name = up.paths.filenames.feat_filt_plot;
        savepath = [up.paths.plots_save_folder, save_name];
        print(h_fig,'-depsc',savepath)
        close all
    end
end


if strfind(up.paths.data_save_folder, 'REST_AND_REC_TEMP')
    %% Create plot of improved precision with no of windows
    fprintf('\n--- Plot of improvement in precision with no of windows ');
    
    %% Load alg names
    load_name = up.paths.filenames.alg_names;
    loadpath = [up.paths.data_save_folder, up.paths.filenames.alg_names, '.mat'];
    load(loadpath, load_name);
    
    %% Load BA data
    load_name = 'BA_results';
    loadpath = [up.paths.data_save_folder, up.paths.filenames.global_BA, '.mat'];
    load(loadpath, load_name);
    rel_BA_res = BA_results.young;
    
    rel_els = find(strcmp(alg_names.sigs, 'ECG'));
    no_wins = alg_names.no_wins(rel_els);
    rel_prec = rel_BA_res.prec.val(rel_els);
    
    % Make Figure
    h_fig = figure('Position', [200, 200, 600, 300]);
    
    % Plot
    ftsize = 12; lwidth = 2;
    f2 = fit(no_wins,rel_prec,'exp2');
    coeffs = coeffvalues(f2);
    fitexp = ( coeffs(1)*exp(coeffs(2)*no_wins) ) + ( coeffs(3)*exp(coeffs(4)*no_wins) );
    plot(no_wins,rel_prec, 'mx', 'LineWidth',lwidth,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10)
    hold on
    plot(no_wins, fitexp, 'r', 'LineWidth', lwidth)
    ylim([8 12])
    xlim([1 max(no_wins)])
    xlabel('Number of windows used to estimate RR', 'FontSize', ftsize)
    ylabel('LOA interval [bpm]', 'FontSize', ftsize)
    set(gca, 'FontSize', ftsize)
    set(gca, 'YTick', 8:12)
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize', [6, 3]);
    set(gcf,'PaperPosition',[0 0 6 3]);
    save_name = up.paths.filenames.temp_prec_plot;
    savepath = [up.paths.plots_save_folder, save_name];
    print(h_fig,'-depsc',savepath)
    close all
end

fprintf('\n--- Making B-A Plots ');

%% Load data from entire study
load_name = up.paths.filenames.win_data;
loadpath = [up.paths.data_save_folder, up.paths.filenames.win_data, '.mat'];
load(loadpath, load_name);

%% Load algorithms names
load_name = up.paths.filenames.alg_names;
loadpath = [up.paths.data_save_folder, up.paths.filenames.alg_names, '.mat'];
load(loadpath, load_name);

%% Load relevant algorithm els
load_name = 'vortal_results';
loadpath = [up.paths.data_save_folder, up.paths.filenames.vortal_results, '.mat'];
load(loadpath, load_name);

%% Identify the most precise algorithm(s) to make plots for
no_plots_per_sig = 1;
for sig = {'ekg', 'ppg'}
    eval([sig{1,1}, '_els = vortal_results.top_ten.' sig{1,1}, '.alg_nos(1:no_plots_per_sig);']);
    eval([sig{1,1}, '_bias = vortal_results.top_ten.' sig{1,1}, '.bias.val(1:no_plots_per_sig);']);
    eval([sig{1,1}, '_lloa = vortal_results.top_ten.' sig{1,1}, '.lloa.val(1:no_plots_per_sig);']);
    eval([sig{1,1}, '_uloa = vortal_results.top_ten.' sig{1,1}, '.uloa.val(1:no_plots_per_sig);']);
end

rel_els = [ekg_els(:); ppg_els(:)];
rel_bias = [ekg_bias(:); ppg_bias(:)];
rel_lloa = [ekg_lloa(:); ppg_lloa(:)];
rel_uloa = [ekg_uloa(:); ppg_uloa(:)];

rel_names = alg_names.names(rel_els);
rel_sigs = alg_names.sigs(rel_els);

%% Generate plots for each algorithm
for rel_el_no = 1 : length(rel_els)
    
    alg_no = rel_els(rel_el_no);
    
    % identify data for this algorithm
    curr_els = win_data.alg_no == alg_no & win_data.snr_log;
    curr_est = win_data.est(curr_els);
    curr_ref = win_data.ref(curr_els);
    curr_bias = rel_bias(rel_el_no);
    curr_lloa = rel_lloa(rel_el_no);
    curr_uloa = rel_uloa(rel_el_no);
    curr_data.ave = nanmean([curr_est(:), curr_ref(:)],2);
    curr_data.error = [curr_est(:) - curr_ref(:)];
    
    % setup figure
    ftsize = 13;
    min_ave_val = 0;            % xlims
    max_ave_val = 36;
    min_error_val = -10;        % ylims
    max_error_val = 10;
    x_edge_int = 0.2;           % resolution of BA plot
    y_edge_int = 0.2;
    edges1 = [min_ave_val : x_edge_int : max_ave_val];       % locations of individual pixels
    edges2 = [min_error_val : y_edge_int : max_error_val];
    
    % Make color plot
    h_fig = figure('Position', [200 200 700 350]);
    A2=hist3([curr_data.ave,curr_data.error],'Edges',{edges1 edges2});
    h = fspecial('gaussian', [10 10], 2);
    y = filter2(h, A2');
    y = y./(max(max(y)));
    imagesc(y),
    c_han = colorbar;
    set(c_han,'YTickMode','manual')
    set(c_han, 'YTick',[0,1], 'YTickLabel', {'Min','Max'})
    
    % mymap = repmat(-1*[-1:0.01:0]', 1,3); colormap(mymap)
    hold on
    
    % X Ticks
    x_tick_int = 10/x_edge_int;
    x_label_int = (edges1(end)-edges1(1))*x_tick_int/length(edges1);
    xlims = xlim;
    edges1_labels = round( edges1(1) : x_label_int : edges1(end) )';
    
    x_spec.axis_lims = xlims;
    x_spec.real_lims = [edges1(1), edges1(end)];
    edges1_ticks = convert_real_to_axis_units(edges1_labels, x_spec);
    
    current_top_of_x_axis = xlims(1) + (0.5*(edges1(2)-edges1(1)));
    % adjust = - current_top_of_x_axis;
    set(gca, 'XTick', edges1_ticks, 'XTickLabel', num2str(edges1_labels))
    
    % Y Ticks
    y_tick_int = 10/y_edge_int;
    y_label_int = (edges2(end)-edges2(1))*y_tick_int/length(edges2);
    ylims = ylim;
    edges2_labels = round( edges2(1) : y_label_int : edges2(end) )';
    
    y_spec.axis_lims = ylims;
    y_spec.real_lims = [edges2(1), edges2(end)];
    edges2_ticks = convert_real_to_axis_units(edges2_labels, y_spec);
    
    current_top_of_y_axis = ylims(1) - (0.5*(edges2(2)-edges2(1)));
    %adjust = - current_top_of_y_axis;
    set(gca, 'YTick', edges2_ticks, 'YTickLabel', num2str(edges2_labels))
    
    % Axis Labels
    
    xlabel('Mean: 0.5(RR^{est} + RR^{ref}),  bpm', 'FontSize', ftsize);
    ylabel('Difference:  (RR^{est} - RR^{ref}),  bpm', 'FontSize', ftsize);
    
    % Bias and LOA Labels
    
    % xlim([min_ave_val, max_ave_val]), ylim([min_error_val, max_error_val])
    %xlims = xlim;
    %bias_adj = BA_data.bias + adjust;
    %loa_adj = [BA_data.lloa, BA_data.uloa] + adjust;
    axis_bias = convert_real_to_axis_units(curr_bias, y_spec);
    axis_lloa = convert_real_to_axis_units(curr_lloa, y_spec);
    axis_uloa = convert_real_to_axis_units(curr_uloa, y_spec);
    plot(xlims, [axis_bias, axis_bias], 'r')
    plot(xlims, [axis_lloa, axis_lloa], '--r')
    plot(xlims, [axis_uloa, axis_uloa], '--r')
    text(mean(xlims),0.5*axis_lloa,['Lower LOA = ' num2str(curr_lloa,2)],'HorizontalAlignment','center','BackgroundColor',[1 1 1], 'FontSize', ftsize);
    text(mean(xlims),1.2*axis_uloa,['Upper LOA = ' num2str(curr_uloa,2)],'HorizontalAlignment','center','BackgroundColor',[1 1 1], 'FontSize', ftsize);
    text(xlims(2),0.8*axis_bias,['Bias = ' num2str(curr_bias,2)],'HorizontalAlignment','right','BackgroundColor',[1 1 1], 'FontSize', ftsize);
    set(gca, 'FontSize', ftsize)
    
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize', [6, 3]);
    set(gcf,'PaperPosition',[0 0 6 3]);
    save_name = ['BA_' rel_sigs{rel_el_no}, rel_names{rel_el_no}];
    savepath = [up.paths.plots_save_folder, save_name];
    print(h_fig,'-depsc',savepath)
    close all
    
end


% %% CDF of abs errors plot
% fprintf('\n--- Making CDF Error Plot');
% 
% % Load data (Best alg)
% load_name = up.paths.filenames.win_data;
% loadpath = [up.paths.data_save_folder, up.paths.filenames.win_data, '.mat'];
% load(loadpath, load_name);
% good_els = win_data.sqi & win_data.snr_log & win_data.alg_no == 203 & ~isnan(win_data.est);
% alg.est = win_data.est(good_els);
% alg.ref = win_data.ref(good_els);
% alg.error = alg.est - alg.ref;
% 
% % Load data (Moderate alg)
% load_name = up.paths.filenames.win_data;
% loadpath = [up.paths.data_save_folder, up.paths.filenames.win_data, '.mat'];
% load(loadpath, load_name);
% good_els = win_data.sqi & win_data.snr_log & win_data.alg_no == 281 & ~isnan(win_data.est);
% mod_alg.est = win_data.est(good_els);
% mod_alg.ref = win_data.ref(good_els);
% mod_alg.error = mod_alg.est - mod_alg.ref;
% 
% % Load data (Imp)
% load_name = up.paths.filenames.win_data;
% loadpath = [up.paths.data_save_folder, up.paths.filenames.win_data_imp, '.mat'];
% load(loadpath, load_name);
% good_els = win_data.sqi & win_data.snr_log & ~isnan(win_data.est) & ~isnan(win_data.est);
% imp.est = win_data.est(good_els);
% imp.ref = win_data.ref(good_els);
% imp.error = imp.est - imp.ref;
% 
% % Normal data
% norm.error = normrnd(0,2.5,1,1000);
% 
% % Make Figure
% ftsize = 12;
% h_fig = figure('Position', [200, 200, 1200, 550]);
% plot(sort(abs(alg.error)),linspace(0,1,length(alg.error))), hold on
% plot(sort(abs(imp.error)),linspace(0,1,length(imp.error)))
% plot(sort(abs(norm.error)),linspace(0,1,length(norm.error)))
% plot(sort(abs(mod_alg.error)),linspace(0,1,length(mod_alg.error)))
% title('Empirical CDF', 'FontSize', ftsize)
% xlabel('Absolute difference (bpm)', 'FontSize', ftsize)
% ylabel('F(x)', 'FontSize', ftsize)
% xlim([0 40])
% set(gca, 'XTick', [0:10, 15:5:40], 'YTick', [0:0.1:0.9, 0.95, 1.0], 'XTickLabel', {'0','','','','','5','','','','','10','15', '20', '25', '30', '35', '40'});
% legend({'a203', 'impedance', 'normal distribution', 'a281'}, 'Location', 'Best')
% grid on
% 
% set(gca, 'FontSize', ftsize)
% save_name = 'Errors CDF';
% savepath = [up.paths.plots_save_folder, save_name];
% savefig(h_fig,savepath)
% close all


fprintf('\n--- Making HR / RR Scatter Plot');

%% Load data from entire study
load_name = up.paths.filenames.win_data;
loadpath = [up.paths.data_save_folder, up.paths.filenames.win_data, '.mat'];
load(loadpath, load_name);
hr = win_data.hr(win_data.comb_log & win_data.ecg_log & win_data.young_log & win_data.alg_no == 1);
rr = win_data.ref(win_data.comb_log & win_data.ecg_log & win_data.young_log & win_data.alg_no == 1);

% Make Figure
h_fig = figure('Position', [200, 200, 900, 400]);

% Plot Histograms
[counts,centers] = hist(rr,20); hold on
counts = counts/max(counts);
scale = 40;
bar(centers,scale*counts, 'c')
[counts,centers] = hist(hr,20);
counts = counts/max(counts);
scale = 5;
barh(centers,scale*counts, 'c')
rr_new = rr; hr_new = hr;

ftsize = 14;
plot(rr_new,hr_new, 'xk')
xlabel('Respiratory Rate [bpm]', 'FontSize', ftsize)
ylabel('Heart Rate [beats per minute]', 'FontSize', ftsize)
set(gca, 'FontSize', ftsize)
set(gca, 'YTick', 0:20:150)
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [9, 4]);
set(gcf,'PaperPosition',[0 0 9 4]);
save_name = up.paths.filenames.hr_rr_scatter;
savepath = [up.paths.plots_save_folder, save_name];
print(h_fig,'-depsc',savepath)
close all

%% Precision and prop of algs using filter / feature

fprintf('\n--- Making Prop of algs, filt/feat plot');

for plot_type = {'est_techs', 'sigs', 'comps'}
    % Load BA data for algs
    if ~strcmp(up.paths.root_data_folder, 'C:\Documents\Data\VORTAL_REST_AND_REC\')
        loadpath = [up.paths.data_save_folder, up.paths.filenames.global_BA, '.mat'];
        load(loadpath);
        rel_res = BA_results.young.prec;
        [~, orders] = sort(rel_res.val);
        rel_res.val = rel_res.val(orders); rel_res.val = rel_res.val(:)';
        rel_res.uci = rel_res.uci(orders); rel_res.uci = rel_res.uci(:)';
        rel_res.lci = rel_res.lci(orders); rel_res.lci = rel_res.lci(:)';
        if strcmp(plot_type{1,1}, 'est_techs')
            rel_els = alg_names.meths.ef(orders)>0;
            rel_res.pos_label = 'Frequency';
            rel_res.neg_label = 'Time';
        elseif strcmp(plot_type{1,1}, 'sigs')
            rel_els = strcmp(alg_names.sigs(orders), 'ECG');
            rel_res.pos_label = 'ECG';
            rel_res.neg_label = 'PPG';
        elseif strcmp(plot_type{1,1}, 'comps')
            rel_res.label1 = 'None';
            rel_res.label2 = 'X_{A1}';
            rel_res.label3 = 'X_{A4}';
            rel_res.label4 = 'E_{F3}';
            rel_els1 = alg_names.meths.xa~= 1 & alg_names.meths.xa~= 4 & alg_names.meths.ef~= 3;
            rel_els2 = alg_names.meths.xa== 1;
            rel_els3 = alg_names.meths.xa== 4;
            rel_els4 = alg_names.meths.ef== 3;
        end
        
    else
        loadpath = 'C:\Users\pc13\Dropbox\VORTAL\VORTAL_theoret_lims_yhvs\2016_Jan_Submission_to_Phys_Meas\Data\Complete_results_table.csv';
        [num, txt, raw] = xlsread(loadpath);
        headers = raw(1,:);
        % Eliminate IP
        rel_col = find(strcmp(headers, 'Signal'));
        bad_row = find(strcmp(raw(:,rel_col),'IP'));
        raw = raw([1:(bad_row-1), (bad_row+1):length(raw(:,1))], :);
        % Eliminate ones which had a problem with the random effects model (failed to converge)
        rel_col = find(strcmp(headers, 'Problems with Random Effects Model'));
        good_rows = ismember(raw(:,rel_col),'.'); good_rows(1) = 1;
        raw = raw(good_rows,:);
        % Extract relevant data
        rel_col = find(strcmp(headers, 'Overall Rank'));
        ranks = cell2mat(raw(2:end, rel_col));
        [~, rel_res.order] = sort(ranks);
        rel_col = find(strcmp(headers, 'Signal'));
        rel_res.sigs = raw(2:end, rel_col);
        rel_col = find(strcmp(headers, 'Algorithm'));
        rel_res.algs = raw(2:end, rel_col);
        rel_col = find(strcmp(headers, '2SD [bpm]'));
        rel_res.val = cell2mat(raw(2:end, rel_col));
        if strcmp(plot_type{1,1}, 'est_techs')
            rel_els = strfind(rel_res.algs, 'Ef');
            rel_els = not(cellfun('isempty', rel_els));
            rel_res.pos_label = 'Frequency';
            rel_res.neg_label = 'Time';
        elseif strcmp(plot_type{1,1}, 'sigs')
            rel_els = strcmp(rel_res.sigs(rel_res.order), 'ECG');
            rel_res.pos_label = 'ECG';
            rel_res.neg_label = 'PPG';
        elseif strcmp(plot_type{1,1}, 'comps')
            rel_res.label1 = 'None';
            rel_res.label2 = 'X_{A1}';
            rel_res.label3 = 'X_{A4}';
            rel_res.label4 = 'E_{F4}';
            temp = strfind(rel_res.algs, 'Xa1'); rel_els2 = ~cellfun(@isempty,temp);
            temp = strfind(rel_res.algs, 'Xa4'); rel_els3 = ~cellfun(@isempty,temp);
            temp = strfind(rel_res.algs, 'Ef4'); rel_els4 = ~cellfun(@isempty,temp);
            rel_els1 = ~rel_els2 & ~rel_els3 & ~rel_els4;
        end
        
    end
    orank = 1 : length(rel_res.val); orank = orank(:)';
    
    % data for hist
    tot = length(orank)+1;
    if strcmp(plot_type{1,1}, 'comps')
        bin_width = (tot/5);
    else
        bin_width = (tot/5);
    end
    bin_ends = 0:bin_width:tot;
    bin_starts = bin_ends(1:(end-1))+1;
    bin_ends = bin_ends(2:end);
    bin_mids = mean([bin_starts(:)'; bin_ends(:)']);
    clear prop*
    if ~strcmp(plot_type{1,1}, 'comps')
        for bin_no = 1 : length(bin_starts)
            prop1(bin_no) = 100*sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no) & rel_els(:))/sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no));
        end
        prop2 = 100 - prop1;
        y_vals = [prop1; prop2];
    else
        for bin_no = 1 : length(bin_starts)
            prop1(bin_no) = 100*sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no) & rel_els1(:))/sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no));
            prop2(bin_no) = 100*sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no) & rel_els2(:))/sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no));
            prop3(bin_no) = 100*sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no) & rel_els3(:))/sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no));
            prop4(bin_no) = 100*sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no) & rel_els4(:))/sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no));
        end
        y_vals = [prop1; prop2; prop3; prop4];
    end
    
    % setup fig
    h_fig = figure('Position', [200 200 900 550]);
    fontsize = 16;
    
    % upper plot
    xlims = [min(orank), max(orank)];
    if xlims(2) == 439
        xlims(2) = 440;
    end
    ylims = [0 50];
    h_prec = subplot(2,1,1, 'Position', [0.2 0.61 0.76 0.38]); hold on
    % fill([orank, fliplr(orank)], [rel_res.lci(:)', fliplr(rel_res.uci(:)')], 0.5*[1,1,1]);
    % hold on,
    plot(orank, rel_res.val, 'k', 'LineWidth', 3),
    ylab = ylabel({'2SD', '[bpm]'}, 'FontSize', fontsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.17, 0.3, 0]);
    xlim(xlims)
    ylim(ylims)
    set(gca, 'XTick', [1, 50:50:tot, tot]);
    xlabel('Overall Rank', 'FontSize', fontsize)
    
    % lower plot (stacked bar)
    h_bar = subplot(2,1,2, 'Position', [0.2 0.10 0.76 0.38]); hold on
    h_bars = bar(bin_mids, y_vals', 'stacked');
    h_bars(1).FaceColor = 0.1*[1,1,1];
    h_bars(2).FaceColor = 0.5*[1,1,1];
    xlim(xlims),
    if ~strcmp(plot_type{1,1}, 'comps')
        ylim([0 100])
    else
        ylim([0 120])
    end
    if strcmp(plot_type{1,1}, 'est_techs')
        ylabel({'% using Frequency- and', 'Time-domain RR estimation'}, 'FontSize', fontsize)
    elseif strcmp(plot_type{1,1}, 'sigs')
        ylabel('% using ECG and PPG', 'FontSize', fontsize)
    elseif strcmp(plot_type{1,1}, 'comps')
        ylab = ylabel({'% using X_{A1},', 'X_{A4}, E_{F4}, or', 'none of these'}, 'FontSize', fontsize, 'Rotation', 0);
        set(ylab, 'Units', 'Normalized', 'Position', [-0.16, 0.33, 0]);
    end
    xlabel('Overall Rank', 'FontSize', fontsize)
    starts = strread(num2str(ceil(bin_starts)),'%s');
    mids = strread(num2str(ceil(bin_mids)),'%s');
    ends = strread(num2str(ceil(bin_ends)),'%s');
    quintiles = strread(num2str(1:length(bin_ends)),'%s');
    for s = 1 : length(ends)
        xtickstr{s} = [starts{s} ' - ' ends{s}];
        %xtickstr{s} = [quintiles{s}, ' (' starts{s} ' - ' ends{s} ')'];
    end
    set(h_bar, 'XTick', bin_mids, 'XTickLabel', xtickstr)
    
    allAxesInFigure = findall(h_fig,'type','axes');
    set(allAxesInFigure, 'Fontsize', fontsize-4, 'xgrid','on')
    
    % annotate
    if ~strcmp(plot_type{1,1}, 'comps')
        text(bin_mids(5),15, rel_res.pos_label, 'HorizontalAlignment', 'Center', 'FontSize', fontsize-4, 'Backgroundcolor', 'w')
        text(bin_mids(5),87, rel_res.neg_label, 'HorizontalAlignment', 'Center', 'FontSize', fontsize-4, 'Backgroundcolor', 'w')
    else
        text(bin_mids(1),50, rel_res.label1, 'HorizontalAlignment', 'Center', 'FontSize', fontsize-4, 'Backgroundcolor', 'w')
        text(bin_mids(5),22, rel_res.label2, 'HorizontalAlignment', 'Center', 'FontSize', fontsize-4, 'Backgroundcolor', 'w')
        text(bin_mids(5),62, rel_res.label3, 'HorizontalAlignment', 'Center', 'FontSize', fontsize-4, 'Backgroundcolor', 'w')
        text(bin_mids(5),95, rel_res.label4, 'HorizontalAlignment', 'Center', 'FontSize', fontsize-4, 'Backgroundcolor', 'w')
    end
    
    % add (a) and (b)
    dim = [.01 .91 .1 .1]; str = '(a)';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',fontsize,'LineStyle','none');
    dim = [.01 .43 .1 .1]; str = '(b)';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',fontsize,'LineStyle','none');
    
    set(gcf,'color','w');
    
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize', [9, 5.5]);
    set(gcf,'PaperPosition',[0 0 9 5.5]);
    save_name = [up.paths.filenames.stacked_bar, '_', plot_type{1,1}];
    savepath = [up.paths.plots_save_folder, save_name];
    print(h_fig,'-depsc',savepath)
    
    close all
    
end

%% New joint stacked bar plot

plot_type = {'comps_est_techs'};
% Load BA data for algs
if strcmp(up.paths.root_data_folder, 'C:\Documents\Data\VORTAL_REST_AND_REC\')
    loadpath = 'C:\Users\pc13\Dropbox\VORTAL\VORTAL_theoret_lims_yhvs\2016_Jan_Submission_to_Phys_Meas\Data\Complete_results_table.csv';
    [num, txt, raw] = xlsread(loadpath);
    headers = raw(1,:);
    % Eliminate IP
    rel_col = find(strcmp(headers, 'Signal'));
    bad_row = find(strcmp(raw(:,rel_col),'IP'));
    raw = raw([1:(bad_row-1), (bad_row+1):length(raw(:,1))], :);
    % Eliminate ones which had a problem with the random effects model (failed to converge)
    rel_col = find(strcmp(headers, 'Problems with Random Effects Model'));
    good_rows = ismember(raw(:,rel_col),'.'); good_rows(1) = 1;
    raw = raw(good_rows,:);
    % Extract relevant data
    rel_col = find(strcmp(headers, 'Overall Rank'));
    ranks = cell2mat(raw(2:end, rel_col));
    [~, rel_res.order] = sort(ranks);
    rel_col = find(strcmp(headers, 'Signal'));
    rel_res.sigs = raw(2:end, rel_col);
    rel_col = find(strcmp(headers, 'Algorithm'));
    rel_res.algs = raw(2:end, rel_col);
    rel_col = find(strcmp(headers, '2SD [bpm]'));
    rel_res.val = cell2mat(raw(2:end, rel_col));
    rel_els = strfind(rel_res.algs, 'Ef');
    rel_els = not(cellfun('isempty', rel_els));
    rel_res.pos_label = 'Frequency';
    rel_res.neg_label = 'Time';
    rel_res.label1 = 'None';
    rel_res.label2 = 'X_{A1}';
    rel_res.label3 = 'X_{A4}';
    rel_res.label4 = 'E_{F4}';
    temp = strfind(rel_res.algs, 'Xa1'); rel_els2 = ~cellfun(@isempty,temp);
    temp = strfind(rel_res.algs, 'Xa4'); rel_els3 = ~cellfun(@isempty,temp);
    temp = strfind(rel_res.algs, 'Ef4'); rel_els4 = ~cellfun(@isempty,temp);
    rel_els1 = ~rel_els2 & ~rel_els3 & ~rel_els4;
    orank = 1 : length(rel_res.val); orank = orank(:)';
    
    % Create first histogram
    tot = length(orank)+1;
    if strcmp(plot_type{1,1}, 'comps')
        bin_width = (tot/5);
    else
        bin_width = (tot/5);
    end
    bin_ends = 0:bin_width:tot;
    bin_starts = bin_ends(1:(end-1))+1;
    bin_ends = bin_ends(2:end);
    bin_mids = mean([bin_starts(:)'; bin_ends(:)']);
    clear prop*
    for bin_no = 1 : length(bin_starts)
        prop1(bin_no) = 100*sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no) & rel_els1(:))/sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no));
        prop2(bin_no) = 100*sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no) & rel_els2(:))/sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no));
        prop3(bin_no) = 100*sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no) & rel_els3(:))/sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no));
        prop4(bin_no) = 100*sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no) & rel_els4(:))/sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no));
    end
    y_vals = [prop1; prop2; prop3; prop4];
    
    % setup fig
    h_fig = figure('Position', [200 200 900 800]);
    fontsize = 16;
    
    % upper plot
    xlims = [min(orank), max(orank)];
    if xlims(2) == 439
        xlims(2) = 440;
    end
    ylims = [0 50];
    h_prec = subplot(3,1,1, 'Position', [0.2 0.74 0.76 0.25]); hold on
    plot(orank, rel_res.val, 'k', 'LineWidth', 3),
    ylab = ylabel({'2SD', '[bpm]'}, 'FontSize', fontsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.17, 0.3, 0]);
    xlim(xlims)
    ylim(ylims)
    set(gca, 'XTick', [1, 50:50:tot, tot]);
    xlabel('Overall Rank', 'FontSize', fontsize)
    
    % mid plot (stacked bar)
    h_bar = subplot(3,1,2, 'Position', [0.2 0.41 0.76 0.25]); hold on
    h_bars = bar(bin_mids, y_vals', 'stacked');
    h_bars(1).FaceColor = 0.1*[1,1,1];
    h_bars(2).FaceColor = 0.5*[1,1,1];
    xlim(xlims),
    ylim([0 120])
    ylab = ylabel({'% using X_{A1},', 'X_{A4}, E_{F4}, or', 'none of these'}, 'FontSize', fontsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.16, 0.33, 0]);
    xlabel('Overall Rank', 'FontSize', fontsize)
    starts = strread(num2str(ceil(bin_starts)),'%s');
    mids = strread(num2str(ceil(bin_mids)),'%s');
    ends = strread(num2str(ceil(bin_ends)),'%s');
    quintiles = strread(num2str(1:length(bin_ends)),'%s');
    for s = 1 : length(ends)
        xtickstr{s} = [starts{s} ' - ' ends{s}];
    end
    set(h_bar, 'XTick', bin_mids, 'XTickLabel', xtickstr)
    % annotate
    text(bin_mids(1),50, rel_res.label1, 'HorizontalAlignment', 'Center', 'FontSize', fontsize-4, 'Backgroundcolor', 'w')
    text(bin_mids(5),25, rel_res.label2, 'HorizontalAlignment', 'Center', 'FontSize', fontsize-4, 'Backgroundcolor', 'w')
    text(bin_mids(5),65, rel_res.label3, 'HorizontalAlignment', 'Center', 'FontSize', fontsize-4, 'Backgroundcolor', 'w')
    text(bin_mids(5),96, rel_res.label4, 'HorizontalAlignment', 'Center', 'FontSize', fontsize-4, 'Backgroundcolor', 'w')
    
    % lower plot (stacked bar)
    for bin_no = 1 : length(bin_starts)
        prop1(bin_no) = 100*sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no) & rel_els(:))/sum(orank(:)>= bin_starts(bin_no) & orank(:)< bin_ends(bin_no));
    end
    prop2 = 100 - prop1;
    y_vals = [prop1; prop2];
    h_bar = subplot(3,1,3, 'Position', [0.2 0.08 0.76 0.25]); hold on
    h_bars = bar(bin_mids, y_vals', 'stacked');
    h_bars(1).FaceColor = 0.1*[1,1,1];
    h_bars(2).FaceColor = 0.5*[1,1,1];
    xlim(xlims),
    ylim([0 100])
    ylab = ylabel({'% using', 'Frequency- or', 'Time-domain', 'RR estimation'}, 'FontSize', fontsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.16, 0.28, 0]);
    xlabel('Overall Rank', 'FontSize', fontsize)
    starts = strread(num2str(ceil(bin_starts)),'%s');
    mids = strread(num2str(ceil(bin_mids)),'%s');
    ends = strread(num2str(ceil(bin_ends)),'%s');
    quintiles = strread(num2str(1:length(bin_ends)),'%s');
    for s = 1 : length(ends)
        xtickstr{s} = [starts{s} ' - ' ends{s}];
    end
    set(h_bar, 'XTick', bin_mids, 'XTickLabel', xtickstr)
    % annotate
    text(bin_mids(5),15, rel_res.pos_label, 'HorizontalAlignment', 'Center', 'FontSize', fontsize-4, 'Backgroundcolor', 'w')
    text(bin_mids(5),87, rel_res.neg_label, 'HorizontalAlignment', 'Center', 'FontSize', fontsize-4, 'Backgroundcolor', 'w')
    
    % add (a) and (b) and (c)
    dim = [.01 .91 .1 .1]; str = '(a)';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',fontsize,'LineStyle','none');
    dim = [.01 .58 .1 .1]; str = '(b)';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',fontsize,'LineStyle','none');
    dim = [.01 .25 .1 .1]; str = '(c)';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',fontsize,'LineStyle','none');
    
    allAxesInFigure = findall(h_fig,'type','axes');
    set(allAxesInFigure, 'Fontsize', fontsize-4, 'xgrid','on')
    
    set(gcf,'color','w');    
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize', [9, 8.0]);
    set(gcf,'PaperPosition',[0 0 9 8.0]);
    save_name = [up.paths.filenames.stacked_bar, '_', plot_type{1,1}];
    savepath = [up.paths.plots_save_folder, save_name];
    print(h_fig,'-depsc',savepath)
    close all
    
end

end

function axis_units = convert_real_to_axis_units(real_units, spec)

axis_units = spec.axis_lims(1) + ( (real_units - spec.real_lims(1)) * (spec.axis_lims(end) - spec.axis_lims(1))/ (spec.real_lims(end)-spec.real_lims(1)));

end