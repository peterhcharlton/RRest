function [respWave] = Wfm(data, up)
%WFM extracts a respiratory signal corresponding to FM using the continuous
%wavelet transform (Morlet wavelet).
%
%	[respWave] = Wfm(data, up)
%
%	Inputs:
%       data            raw signal data
%       up              universal parameters structure
%
%	Outputs:
%       respWave        a respiratory signal
%

%% Setup
data.t = (1/data.fs)*(1:length(data.v));
data.v = data.v;

%% Downsample to make processing possible
up.paramSet.filt_resample_fs = 25;   % temporarily change downsample freq (perhaps to 50) to increase resolution?
d_s = downsample_data(data, up);
up.paramSet.filt_resample_fs = 25;

%% CWT
% Specify characteristics
s0  = 6/d_s.fs;  % smallest scale
ds = 0.001; % spacing between scales
NbSc = 3000; % number of scales
SCA = {s0,ds,NbSc, 'lin'}; % specify scales
cwtstruct = cwtft({d_s.v, 1/d_s.fs},'scales',SCA);
scales = cwtstruct.scales;
F = scal2frq(scales,cwtstruct.wav,1/d_s.fs); F = d_s.fs./F;
%contour(d_s.t,F,real(cwtstruct.cfs));
%xlabel('Seconds'); ylabel('Pseudo-frequency');
%hold on
mag_mat = cwtstruct.cfs;
rel_rows = F >= up.paramSet.hr_range(1)/60 & F <= up.paramSet.hr_range(2)/60;
rel_mag_mat = mag_mat(rel_rows, :);
[rel_mags, rel_els] = max(rel_mag_mat);
%plot(d_s.t,F(rel_els), 'r', 'LineWidth', 3)
[fm_sig, am_sig] = deal(d_s);
fm_sig.v = F(rel_els);
am_sig.v = abs(rel_mags);

%% Downsample result to usual downsample freq
fm_sig = downsample_data(fm_sig, up);
am_sig = downsample_data(am_sig, up);

%% Store
respWave = fm_sig;

end