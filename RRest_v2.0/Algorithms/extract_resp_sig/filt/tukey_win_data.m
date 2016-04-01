function d_s_win = tukey_win_data(d_s, up)
%TUKEY_WIN_DATA extracts respiratory activity using a BPF (PC's implementation)
%
%	d_s_win = tukey_win_data(data, up)
%
%	Inputs:
%       data            raw signal data
%       up              universal parameters structure
%
%	Outputs:
%       d_s_win             a signal after Tukey windowing
%

%% Window Signal to reduce edge effects

duration_of_signal = d_s.t(end) - d_s.t(1);
prop_of_win_in_outer_regions = 2*up.paramSet.tukey_win_duration_taper/duration_of_signal;
tukey_win = tukeywin(length(d_s.v), prop_of_win_in_outer_regions);
d_s_win = d_s;    % copy time and fs
d_s_win.v = detrend(d_s.v(:)).*tukey_win(:);

end
