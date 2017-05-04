function [respWave] = BFi(data, up)
%BFi extracts respiratory activity using a BPF
%
%	[respWave] = BFi(data, up)
%
%	Inputs:
%       data            raw signal data
%       up              universal parameters structure
%
%	Outputs:
%       respWave        a respiratory signal after bpf
%

%% Setup
data.t = (1/data.fs)*(1:length(data.v));
data.v = data.v;

%% Downsample
respWave = downsample_data(data, up);

end
