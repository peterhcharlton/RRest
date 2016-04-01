function [respWave] = CCF(data, up)
%CCF extracts respiratory signal using the centered-correntropy function.
%
%	[respWave] = CCF(data, up)
%
%	Inputs:
%       data            raw signal data
%       up              universal parameters structure
%
%	Outputs:
%       respWave        a respiratory signal after ccf
%

%% Setup
data.t = (1/data.fs)*(1:length(data.v));
data.v = data.v;

%% Downsample
d_s = downsample_data(data, up);

%% Compute Centred Correntropy Function on the pre-processed signal
V_c = corrsd(d_s.v);
respWave = d_s;
respWave.v = V_c;
respWave.t = d_s.t(1:length(respWave.v));

end

function V_c = corrsd(s)
%% Compute Centred Correntropy Function
% N, no of samples in window
N = length(s);
% Calculate the length scale for the Gaussian Kernel, sigma
sigma = calc_sigma(s);

% Find distances and kernel density estimates, kdes
distances = nan(N,N);
for n = 2 : N
    rel_distances_els = 1 : (n-1);
    rel_signal_els = n - rel_distances_els;
    
    distances(rel_distances_els, n) = s(n) - s(rel_signal_els);
end

% kernel density estimates, kdes
kdes = nan(N,N);
for n = 2 : N
    rel_kdes_els = 1 : (n-1);
    
    kdes(rel_kdes_els, n) = gaussian_kernel(distances(rel_kdes_els,n), sigma);
end

% Find summed and mean kdes for each value of m
summed_kdes = nan(N-1,1);
mean_kdes = nan(N-1,1);
for m = 1 : (N-1)
    rel_kdes = (m+1):N;   % This is simply those which were calculated
    summed_kdes(m) = sum(kdes(m,rel_kdes));
    mean_kdes(m) = mean(kdes(m,rel_kdes));
end

% Now we have the correntropy function
V = mean_kdes;

% Correntropy mean:
V_bar = (1/(N^2))*sum(V);

% Centred correntropy function:
V_c = V - V_bar;

end

function sigma = calc_sigma(s)

% N, no of samples in window
N = length(s);

% Calculate sigma
% uses Silverman's rule of density estimation, as described in:
% A. Garde (2010), �Correntropy-based spectral characterization of
% respiratory patterns in patients with chronic heart failure.,�
% IEEE Trans. Biomed. Eng., vol. 57, no. 8, pp. 1964�72, Aug. 2010.

% used to be:
% [q1 q2 q3 fu fl ou ol] = quartile(s); % find quartiles with which to calculate iqr
% interqr = q3-q1;
% now:
interqr = iqr(s);
A = min([std(s), 1.34*interqr]);
sigma = 0.9*A*(N^-0.2);

end

function kappa = gaussian_kernel(distances, sigma)

%% Calculate Kappa using Gaussian kernel
kappa = (1/(sqrt(2*pi)*sigma)).*exp( -1*(distances.^2)./ (2*(sigma^2)) );

end