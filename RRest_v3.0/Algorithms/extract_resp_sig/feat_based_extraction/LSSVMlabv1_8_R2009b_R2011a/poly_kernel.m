function x = poly_kernel(a,b, d)
% polynomial kernel function for implicit higher dimension mapping
%
%  X = poly_kernel(a,b,[t,degree])
%
% 'a' can only contain one datapoint in a row, 'b' can contain N
% datapoints of the same dimension as 'a'. 
%
% x = (a*b'+t^2).^degree;
% 
% see also:
%    RBF_kernel, lin_kernel, MLP_kernel, trainlssvm, simlssvm
%

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

if length(d)>1, 
    t=d(1);
    d=d(2); 
else
    d = d(1);t=1; 
end
d = (abs(d)>=1)*abs(d)+(abs(d)<1); % >=1 !!

x=(b*a'+ t^2).^d;



