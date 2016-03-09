function identify_subj_wins( subj, up )
%IDENTIFY_SUBJ_WINS Identifies windows in the respiratory signals on which
%to perform RR estimation. 
%
%   inputs -    subj        - the number of the subject (according to up.paramSet.subj_list)
%               up          .paramSet.winLeng - the duration of each window in secs
%               up          .paramSet.buffer_period - the number of secs to ignore at the start of each resp sig (since the filters might not have stabilised).
%               up          .paramSet.winStep - the number of secs between each consecutive window
%   outputs -   wins        .t - vector of start times
%

%% Load data if it already exists
save_name = 'wins';
savepath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.win_timings, '.mat'];
if exist(savepath, 'file')
    return
end

%% Load respiratory signals for this subject
respSigs = load([up.paths.data_save_folder, num2str(subj), up.paths.filenames.respSigs]);

%% Identify start and end times of each respiratory signal
respSig_names = fieldnames(respSigs);
[timings.start, timings.end] = deal(nan(length(respSig_names),1));
for sig_no = 1 : length(respSig_names)
    eval(['rel_data = respSigs.' respSig_names{sig_no} '.t;']);
    timings.start(sig_no) = rel_data(1);
    timings.end(sig_no) = rel_data(end);
end

latest_start = max(timings.start);
earliest_end = min(timings.end);

%% Define windows
duration_of_one_win = up.paramSet.winLeng;
no_of_secs_bet_con_wins = up.paramSet.winStep;
gap_between_win_starts = duration_of_one_win - no_of_secs_bet_con_wins;
wins.t_start = latest_start : gap_between_win_starts : (earliest_end - duration_of_one_win);
wins.t_end = wins.t_start + duration_of_one_win;

%% Save
save_or_append_data

end