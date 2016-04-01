function d_s = downsample_data(data, up)
%DOWNSAMPLE_DATA extracts respiratory activity using a BPF (PC's implementation)
%
%	d_s = downsample_data(data, up)
%
%	Inputs:
%       data            raw signal data
%       up              universal parameters structure
%
%	Outputs:
%       d_s             signal after downsampling
%

%% Downsample
downsample_freq = up.paramSet.filt_resample_fs;
d_s.fs = downsample_freq;
d_s.t = downsample(data.t, data.fs/downsample_freq);
d_s.v = decimate(data.v, data.fs/downsample_freq);

end
