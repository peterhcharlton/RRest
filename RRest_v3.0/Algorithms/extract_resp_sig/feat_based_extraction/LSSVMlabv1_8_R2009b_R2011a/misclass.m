function [perc,n,which] = misclass(Y,Yest)
% The rate of misclassifications.
%
% '[perc,n,which] = misclass(Y,Yest)'
%
%  'Y' contains the real class labels; 
%  'Yest' contains the estimated class labels;
%  'perc' is the rate of misclassifications (between 0 and 1); 
%  'n' is the number of misclassifications;
%  'which' contains the indices of the misclassificated instances
%     (the first column gives the row, the second the column index)
%
%
% see also:
%    validate, mse, linf, medae, mae

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

n = sum(sum(Y~=Yest));
perc = n/numel(Y);
[I,J] = find(Y~=Yest);
which = [J I];