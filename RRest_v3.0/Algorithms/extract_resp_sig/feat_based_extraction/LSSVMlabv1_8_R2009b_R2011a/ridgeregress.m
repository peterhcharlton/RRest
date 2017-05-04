function [w,b,Yt] = ridgeregress(X,Y,gam, Xt)
% Linear ridge regression
% 
% >> [w, b]     = ridgeregress(X, Y, gam)
% >> [w, b, Yt] = ridgeregress(X, Y, gam, Xt)
% 
% Ordinary Least squares with a regularization parameter (gam).
% 
% Full syntax
% 
% >> [w, b, Yt] = ridgeregress(X, Y, gam, Xt)
% 
% Outputs    
%   w     : d x 1 vector with the regression coefficients
%   b     : bias term
%   Yt(*) : Nt x 1 vector with predicted outputs of test data
% Inputs    
%   X     : N x d matrix with the inputs of the training data
%   Y     : N x 1 vector with the outputs of the training data
%   gam   : Regularization parameter
%   Xt(*) : Nt x d matrix with the inputs of the test data
% 
% See also:
% bay_rr,bay_lssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

  
if size(X,1)~=size(Y,1),
  error('X and Y need to have the same number of data points');
end
if size(Y,2)~=1,
  error('Only handling one-dimensional output');
end
if nargin==4 & size(Xt,2)~=size(X,2),
  error('Training input and test inputs need to have the same dimension');
end
  
[nD,nx] = size(X);
if nx>nD, warning('dim datapoints larger than number of datapoints...');end


Xe = [X ones(nD,1)];
%H = [ Xe'*Xe + gam^-1*[eye(nx) zeros(nx,1); zeros(1,nx+1)]];
H = Xe'*Xe + inv(gam).*eye(nx+1);

sol = pinv(H)*Xe'*Y;
w = sol(1:end-1);
b = sol(end);


if nargin<4, return; end
Yt = Xt*w+b;
  
  
   
