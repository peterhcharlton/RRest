function x = lin_kernel(a,b)
% kernel function for implicit higher dimension mapping, based on
% the standard inner-product
%
%   x = lin_kernel(a,b)
%
% 'a' can only contain one datapoint in a row, 'b' can contain N
% datapoints of the same dimension as 'a'. 
%
% see also:
%    poly_kernel, RBF_kernel, MLP_kernel, trainlssvm, simlssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


x = zeros(size(b,1),1);
for i=1:size(b,1),
  x(i,1) = a*b(i,:)';
end