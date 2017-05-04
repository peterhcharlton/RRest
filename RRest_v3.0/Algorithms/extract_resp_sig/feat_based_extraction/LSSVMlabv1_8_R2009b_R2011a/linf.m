function [perf,which] = linf(e)
% L infinity norm of the residuals
% 
%  >> perf = linf(E);
%
% see also:
%    mse, mae, medae, trimmedmse

  
% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

[m1,w1] = max(abs(e));
[perf,w2] = max(m1);
which = [w1(w2) w2];