function feature_measurement_for_breathing_rate_estimation
% feature_measurement_for_breathing_rate_estimation creates a plot of
% feature measurements of the ECG and PPG signals for breathing rate
% estimation.
%
%               feature_measurement_for_breathing_rate_estimation
%
%	This file creates the image originally published in:
%           Charlton P.H. et al. Breathing Rate Estimation from the
%           Electrocardiogram and Photoplethysmogram: A Review,
%           IEEE Reviews in Biomedical Engineering, 2017
%           DOI: https://doi.org/10.1109/RBME.2017.2763681
%   Please cite this publication when using this image.
%
%   The image is adapted from Fig. 26.6 in:
%           Charlton P.H. et al. Waveform analysis to estimate respiratory
%           rate, in Secondary Analysis of Electronic Health Records.
%           Springer International Publishing, pp. 37790, 2016. 
%           DOI: https://doi.org/10.1007/978-3-319-43742-2_26
%   
%   Output:
%       PNG and PDF images in the same folder as this script
%
%   Pre-requisites:
%       This script uses the "arrow" function available at:
%           https://uk.mathworks.com/matlabcentral/fileexchange/278-arrow
%       The location of this function must be in the Matlab path for this
%       script to work.
%     
%   Comments, Questions, Criticisms, Feedback, Contributions:
%       See: http://peterhcharlton.github.io/RRest/contributions.html
%
%   Licence:
%       please see the accompanying file named "LICENSE"
%

%% Check whether the "arrow" function is available
if ~exist('arrow.m', 'file')
    error('This script requires the "arrow" function, available at: \n%s', 'https://uk.mathworks.com/matlabcentral/fileexchange/278-arrow')
end

%% setup path
% find out where this script is:
path_of_this_script = mfilename('fullpath');
[location_of_this_script,~,~] = fileparts(path_of_this_script);
cd(location_of_this_script)

%% Settings
no_beats = 511;
sig_range = 1;
f_r = 0.25; w_r = 2*pi*f_r;
f_c = 1.2; w_c = 2*pi*f_c;
fs = 100;
a_r = 0.22;

%% Load data

ecg_data.t = 0:0.01:1;
ecg_data.v = [0.0628929206915458,0.0628929206915458,0.0587000506488106,0.0607963906957733,0.0587000506488106,0.0628929206915458,0.0587000506488106,0.0649893007277317,0.0587000506488106,0.0649893906993127,0.0628929206915458,0.0628929206915458,0.0587001706021793,0.0607964806721212,0.0754714208624968,0.0670857907342811,0.0712786607770163,0.0607965306526915,0.0587000506488106,0.0712790306332367,0.0796641009790649,0.0880504708590479,0.100628751076163,0.0985321711111414,0.0922430109906925,0.0817607109175321,0.0796644008624868,0.0859536009673920,0.0838574908197314,0.0943395309891593,0.0964358810334277,0.0964358810334277,0.0796640410023805,0.138365121315293,0.553456376600879,1,0.964354102164553,0.276729142911456,0,0.0251570603185868,0.0628930206526865,0.0649893907034836,0.0712784908430773,0.0733751507770112,0.0712787207537007,0.0796644008624868,0.0796644008624868,0.0880501409479572,0.0838572908974501,0.0943395809596000,0.0964358810334277,0.106917921192726,0.109014511153862,0.115303891188820,0.121593101289839,0.129978841375310,0.138364581460780,0.153039681588981,0.167714841693866,0.186582621942805,0.218028912350468,0.247379232560238,0.285115292855479,0.322850853343825,0.352200673749083,0.375261713885036,0.389936784024894,0.394129654067629,0.373165483784006,0.329140363329457,0.259957612777821,0.190775711926268,0.138365241204309,0.0922432209090877,0.0754715108275234,0.0545071806060754,0.0545071806060754,0.0461214405206049,0.0545071806060754,0.0566034506718976,0.0628929206915458,0.0628929206915458,0.0712786607770163,0.0712786607770163,0.0754715308197515,0.0754715308197515,0.0796644008624868,0.0796644008624868,0.0796644008624868,0.0796644008624868,0.0754716807614624,0.0712786607770163,0.0712786607770163,0.0754715308197515,0.0670859806604482,0.0628929206915458,0.0670857907342811,0.0607963906957734,0.0670857907342811,0.0628929206915458,0.0712786607770163];

ppg.t = 0:0.01:1;
ppg.v = [0,0.00927542939852481,0.0372488322480372,0.0896445546165244,0.169004189004819,0.272647048523832,0.393762147326141,0.522824548449564,0.648130256209061,0.758626219281831,0.847819546835619,0.914400499067893,0.959234915642351,0.985437381517778,0.997714838452157,1,0.995011418768304,0.985125871133712,0.971761461825674,0.955456925699569,0.937412807053363,0.917629320396664,0.896682390477446,0.875479725915258,0.853916479121250,0.832170367477851,0.810423915529333,0.788209093110292,0.764908967193013,0.739300355899694,0.712381664572962,0.683143829496475,0.651587233363082,0.618067496536693,0.583173995134062,0.547292623927837,0.511984733125585,0.478482376769321,0.447658087846185,0.420523059079013,0.397447373925289,0.378897500664425,0.366497617029747,0.359654515162849,0.357484309046336,0.359638663309276,0.364692117662492,0.372113594720006,0.380055318961484,0.388913036314599,0.399242502208559,0.406355757469157,0.408744202083830,0.407861780008873,0.404510362335276,0.398839763639832,0.390849762384774,0.382094285048579,0.372173831593098,0.361197192580294,0.349640925944862,0.337504695928278,0.325000637378979,0.312496521432379,0.300276611142535,0.287995625629235,0.275491584041999,0.263380083238680,0.251419743047758,0.239459173263126,0.227498880486464,0.215538282837730,0.203643944307285,0.192226987553774,0.180948415515039,0.170075343432093,0.159202176822953,0.148329061636135,0.137738867025047,0.127409366648891,0.117434993958546,0.108040650505145,0.0992259833676540,0.0909912956489715,0.0833362360009967,0.0775118991588359,0.0715432265228746,0.0656005086235839,0.0601638844825797,0.0549468564227474,0.0507434502364092,0.0467557643993294,0.0426247533240844,0.0390735410794610,0.0355210529792382,0.0310620642768787,0.0254434048267722,0.0182660054583478,0.0102196996141933,0.00347880408741879,0.00152249485121779];

sin_wave.t = ppg.t;
sin_wave.v = sin(w_c*ppg.t);

%% Find plain data
for sig = {'ppg', 'ecg_data', 'sin_wave'}
    % Extract data
    eval(['data = ' sig{1,1} ';']);
    % Normalise data
    data.v = sig_range * (data.v - min(data.v))/(range(data.v));
    % Multiply data
    data.seg.v = [ data.v, repmat(data.v(2:end), [1,no_beats-1]) ];
    data.seg.t = 0:(1/fs):no_beats;
    % Re-insert data
    eval([sig{1,1} ' = data;']);
end

%% Find modulated data
for sig = {'ppg', 'ecg_data', 'sin_wave'}
    % Extract data
    eval(['data = ' sig{1,1} ';']);
    % BW:
    data.bw.v = data.seg.v + a_r*sin(w_r*data.seg.t);
    data.bw.t = data.seg.t;
    % AM:
    data.am.t = data.seg.t;
    % data.am.v = data.seg.v.*(1+a_r*sin(w_r*data.seg.t));
    data.am.v = detrend(data.seg.v).*((1/a_r)+cos(w_r*data.seg.t));
    % FM:
    temp.t = data.seg.t + a_r*sin(w_r*data.seg.t);
    temp.v = data.seg.v;
    data.fm.t = data.seg.t;
    data.fm.v = interp1(temp.t, temp.v, data.fm.t);
    % Re-insert data
    eval([sig{1,1} ' = data;']);
end


%% Make Figure demonstrating feature measurement

% setup figure
duration = 6;
start_t = 4.8;
paper_size = [200, 200, 900, 600];
h_fig = figure('Position',paper_size);
fontsize = 16;
width = 0.4;
height = 0.28;
marksize = 6;
lwidth = 2;
ylabpos = [-0.09, 0.4];

% Plot ECG
ecg_data.pk_els = 36:100:length(ecg_data.seg.v);
ecg_data.tr_els = ecg_data.pk_els-3;
temp = ecg_data;
ecg_col = [0.35, 0.45, 0.2];

subpos = [0.08, 0.70, width, height];
ax(1) = subplot('Position', subpos);
data = temp.bw; rel_els = find(data.t>= start_t & data.t<=start_t+duration);
plot(data.t(rel_els), data.v(rel_els), 'Color', ecg_col, 'LineWidth', lwidth),
hold on, plot(data.t(ecg_data.pk_els), data.v(ecg_data.pk_els), 'ok', 'LineWidth', lwidth, 'MarkerSize', marksize, 'MarkerEdgeColor','k','MarkerFaceColor','k')
plot(data.t(ecg_data.tr_els), data.v(ecg_data.tr_els), 'ok', 'LineWidth', lwidth, 'MarkerSize', marksize, 'MarkerEdgeColor','k','MarkerFaceColor','k')
temp2.v = mean([data.v(ecg_data.tr_els); data.v(ecg_data.pk_els)],1);
temp2.t = mean([data.t(ecg_data.tr_els); data.t(ecg_data.pk_els)],1);
plot(temp2.t, temp2.v, 'o-r', 'LineWidth', lwidth, 'MarkerSize', marksize, 'MarkerEdgeColor','r','MarkerFaceColor','r')
xlim([min(data.t(rel_els)), max(data.t(rel_els))]),
ylim([min(data.v(rel_els))-0.2, max(data.v(rel_els))+0.2]),
set(gca, 'XTick', [], 'YTick', []),
ylab = ylabel({'ECG ', '(BW) '}, 'Rotation', 0, 'Color', 'k', 'FontSize', fontsize);
set(ylab, 'Units', 'Normalized', 'Position', ylabpos);

subpos = [0.08, 0.38, width, height];
ax(2) = subplot('Position', subpos);
data = temp.am; rel_els = find(data.t>= start_t & data.t<=start_t+duration);
plot(data.t(rel_els), data.v(rel_els), 'Color', ecg_col, 'LineWidth', 2),
hold on, plot(data.t(ecg_data.pk_els), data.v(ecg_data.pk_els), 'ok', 'LineWidth', lwidth, 'MarkerSize', marksize, 'MarkerEdgeColor','k','MarkerFaceColor','k')
plot(data.t(ecg_data.tr_els), data.v(ecg_data.tr_els), 'ok', 'LineWidth', lwidth, 'MarkerSize', marksize, 'MarkerEdgeColor','k','MarkerFaceColor','k')
temp2.t = mean([data.t(ecg_data.tr_els); data.t(ecg_data.pk_els)],1)-0.15;
temp2.top = data.v(ecg_data.pk_els);
temp2.len = data.v(ecg_data.pk_els) - data.v(ecg_data.tr_els);
xlim([min(data.t(rel_els)), max(data.t(rel_els))]),
ylim([min(data.v(rel_els))-0.7, max(data.v(rel_els))+0.9]),
set(gca, 'XTick', [], 'YTick', []),
ylab = ylabel({'ECG ', '(AM) '}, 'Rotation', 0, 'Color', 'k', 'FontSize', fontsize);
set(ylab, 'Units', 'Normalized', 'Position', ylabpos);
for s = 1:length(ecg_data.pk_els)
    if temp2.t(s) >= start_t && temp2.t(s) <= start_t+duration
        arrow([temp2.t(s), temp2.top(s)-temp2.len(s)], [temp2.t(s), temp2.top(s)], 'Width', 3, 'Ends', 'both', 'FaceColor', [1,0,0], 'EdgeColor', 'none');
    end
end

subpos = [0.08, 0.06, width, height];
ax(3) = subplot('Position', subpos);
data = temp.fm; rel_els = find(data.t>= start_t & data.t<=start_t+duration);
ecg_data.pk_els = [318,448,555,625,718,848,955,1025,1118];
ecg_data.tr_els = ecg_data.pk_els-3;
plot(data.t(rel_els), data.v(rel_els), 'Color', ecg_col, 'LineWidth', 2),
hold on, plot(data.t(ecg_data.pk_els), data.v(ecg_data.pk_els), 'ok', 'LineWidth', lwidth, 'MarkerSize', marksize, 'MarkerEdgeColor','k','MarkerFaceColor','k')
xlim([min(data.t(rel_els)), max(data.t(rel_els))]),
ylim([min(data.v(rel_els))-0.2, max(data.v(rel_els))+0.2]),
set(gca, 'XTick', [], 'YTick', []),
xlabel('Time', 'FontSize', fontsize),
ylab = ylabel({'ECG ', '(FM) '}, 'Rotation', 0, 'Color', 'k', 'FontSize', fontsize);
temp2.t = data.t(ecg_data.pk_els);
temp2.top = data.v(ecg_data.pk_els);
set(ylab, 'Units', 'Normalized', 'Position', ylabpos);
for s = 1:length(ecg_data.pk_els)
    if s>1 && temp2.t(s-1) >= start_t && temp2.t(s) <= start_t+duration
        arrow([temp2.t(s-1), mean(temp2.top)], [temp2.t(s), mean(temp2.top)], 'Width', 3, 'Ends', 'both', 'FaceColor', [1,0,0], 'EdgeColor', 'none');
    end
end

% Plot PPG
ppg.pk_els = 16:100:length(ppg.seg.v);
ppg.tr_els = ppg.pk_els-15;
temp = ppg;
ppg_col = [0.3, 0.3, 1];

subpos = [0.58, 0.70, width, height];
ax(1) = subplot('Position', subpos);
data = temp.bw; rel_els = find(data.t>= start_t & data.t<=start_t+duration);
plot(data.t(rel_els), data.v(rel_els), 'Color', ppg_col, 'LineWidth', lwidth),
hold on, plot(data.t(ppg.pk_els), data.v(ppg.pk_els), 'ok', 'LineWidth', lwidth, 'MarkerSize', marksize, 'MarkerEdgeColor','k','MarkerFaceColor','k')
plot(data.t(ppg.tr_els), data.v(ppg.tr_els), 'ok', 'LineWidth', lwidth, 'MarkerSize', marksize, 'MarkerEdgeColor','k','MarkerFaceColor','k')
temp2.v = mean([data.v(ppg.tr_els); data.v(ppg.pk_els)],1);
temp2.t = mean([data.t(ppg.tr_els); data.t(ppg.pk_els)],1);
plot(temp2.t, temp2.v, 'o-r', 'LineWidth', lwidth, 'MarkerSize', marksize, 'MarkerEdgeColor','r','MarkerFaceColor','r')
xlim([min(data.t(rel_els)), max(data.t(rel_els))]),
ylim([min(data.v(rel_els))-0.2, max(data.v(rel_els))+0.2]),
set(gca, 'XTick', [], 'YTick', []),
ylab = ylabel({'PPG ', '(BW) '}, 'Rotation', 0, 'Color', 'k', 'FontSize', fontsize);
set(ylab, 'Units', 'Normalized', 'Position', ylabpos);

subpos = [0.58, 0.38, width, height];
ax(2) = subplot('Position', subpos);
data = temp.am; rel_els = find(data.t>= start_t & data.t<=start_t+duration);
plot(data.t(rel_els), data.v(rel_els), 'Color', ppg_col, 'LineWidth', 2),
hold on, plot(data.t(ppg.pk_els), data.v(ppg.pk_els), 'ok', 'LineWidth', lwidth, 'MarkerSize', marksize, 'MarkerEdgeColor','k','MarkerFaceColor','k')
plot(data.t(ppg.tr_els), data.v(ppg.tr_els), 'ok', 'LineWidth', lwidth, 'MarkerSize', marksize, 'MarkerEdgeColor','k','MarkerFaceColor','k')
temp2.t = mean([data.t(ppg.tr_els); data.t(ppg.pk_els)],1)-0.15;
temp2.top = data.v(ppg.pk_els);
temp2.len = data.v(ppg.pk_els) - data.v(ppg.tr_els);
xlim([min(data.t(rel_els)), max(data.t(rel_els))]),
ylim([min(data.v(rel_els))-0.7, max(data.v(rel_els))+0.9]),
set(gca, 'XTick', [], 'YTick', []),
ylab = ylabel({'PPG ', '(AM)'}, 'Rotation', 0, 'Color', 'k', 'FontSize', fontsize);
set(ylab, 'Units', 'Normalized', 'Position', ylabpos);
for s = 1:length(ppg.pk_els)
    if temp2.t(s) >= start_t && temp2.t(s) <= start_t+duration
        arrow([temp2.t(s), temp2.top(s)-temp2.len(s)], [temp2.t(s), temp2.top(s)], 'Width', 3, 'Ends', 'both', 'FaceColor', [1,0,0], 'EdgeColor', 'none');
    end
end

subpos = [0.58, 0.06, width, height];
ax(3) = subplot('Position', subpos);
ppg.pk_els = [294,421,537,611,694,821, 937,1011,1094,1221];
ppg.tr_els = ppg.pk_els-3;
data = temp.fm; rel_els = find(data.t>= start_t & data.t<=start_t+duration);
plot(data.t(rel_els), data.v(rel_els), 'Color', ppg_col, 'LineWidth', 2),
hold on, plot(data.t(ppg.pk_els), data.v(ppg.pk_els), 'ok', 'LineWidth', lwidth, 'MarkerSize', marksize, 'MarkerEdgeColor','k','MarkerFaceColor','k')
xlim([min(data.t(rel_els)), max(data.t(rel_els))]),
ylim([min(data.v(rel_els))-0.2, max(data.v(rel_els))+0.2]),
set(gca, 'XTick', [], 'YTick', []),
xlabel('Time', 'FontSize', fontsize),
ylab = ylabel({'PPG ', '(FM) '}, 'Rotation', 0, 'Color', 'k', 'FontSize', fontsize);
temp2.t = data.t(ppg.pk_els);
temp2.top = data.v(ppg.pk_els);
set(ylab, 'Units', 'Normalized', 'Position', ylabpos);
for s = 1:length(ppg.pk_els)
    if s>1 && temp2.t(s-1) >= start_t && temp2.t(s) <= start_t+duration
        arrow([temp2.t(s-1), mean(temp2.top)], [temp2.t(s), mean(temp2.top)], 'Width', 3, 'Ends', 'both', 'FaceColor', [1,0,0], 'EdgeColor', 'none');
    end
end

%% Print to file
savepath = [location_of_this_script, filesep, 'resp_mod_measures'];

% PNG
print(h_fig,'-dpng',savepath)

% PDF
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [paper_size(3)/100, paper_size(4)/100]);
set(gcf,'PaperPosition',[0 0 paper_size(3)/100 paper_size(4)/100]);
print(h_fig,'-dpdf',savepath)

% EPS
% I found that using the save command on the figure window worked best.

close all

end
