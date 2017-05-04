function perf=mse(e)
%
% calculate the mean squared error of the given errors
% 
%  'perf = mse(E);'
%
% see also:
%    mae, linf, trimmedmse
%

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


perf = sum(sum(e.^2)) / numel(e);