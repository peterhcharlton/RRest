function perf=medae(e)
% calculate the median absolute error of the given errors
% 
%  >> perf = medae(E);
%
% see also:
%    mse, mae, linf, trimmedmse
%

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

perf = median(abs(reshape(e,numel(e),1)));
