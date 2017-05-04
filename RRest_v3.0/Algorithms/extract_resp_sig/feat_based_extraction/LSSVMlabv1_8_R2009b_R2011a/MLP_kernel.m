function x = MLP_kernel(a,b, par)
% Multi Layer Perceptron kernel function for implicit higher dimension mapping
%
%   x = MLP_kernel(a,b,[s,t])
% 
% 'a' can only contain one datapoint in a row, 'b' can contain N
% datapoints of the same dimension as 'a'. 
%
%   x = tanh(s*a'b+t^2)
%
% see also:
%    poly_kernel, lin_kernel, RBF_kernel, trainlssvm, simlssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

if length(par)==1, par(2) = 1; end 
x = zeros(size(b,1),1);
for i=1:size(b,1),
  dp = a*b(i,:)';
  x(i,1) = tanh(par(1)*dp + par(2)^2);
end