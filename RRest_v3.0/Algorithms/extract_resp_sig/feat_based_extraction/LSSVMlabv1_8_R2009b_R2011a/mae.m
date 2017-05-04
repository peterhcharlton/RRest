function perf=mae(e)
%
% calculate the absolute error of the given errors
% 
%  'perf = mae(E);'
%
% see also:
%    mse, linf
%

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


perf = sum(sum(abs(e))) / numel(e);