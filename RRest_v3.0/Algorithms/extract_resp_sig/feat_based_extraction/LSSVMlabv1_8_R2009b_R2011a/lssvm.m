function [yp,alpha,b,gam,sig2,model] = lssvm(x,y,type,varargin)
%
% one line LS-SVM calculation
% 
% >> [yp,alpha,b,gam,sig2,model] = lssvm(x,y,type,varargin)
%
% x is the input data, N x d, (can be uni- or multivariate) and y, N x 1, is the response 
% variable 
%
% syntax:  
%       RBF-kernel is used with standard simplex method
%       yp = lssvm(x,y,'f') 
%  
%       lin/poly/RBF is used with standard simplex      
%       yp = lssvm(x,y,'f',kernel)  
%
%  output:
%       yp    : N x 1 vector of predicted outputs
%       alpha : N x 1 vector of lagrange multipliers of the LS-SVM
%       b     : LS-SVM bias term
%       gam   : tuned regularization constant
%       sig2  : squared tuned kernel bandwidth
%       model : object oriented interface of the LS-SVM
%
% See also:
%   trainlssvm, simlssvm, crossvalidate, leaveoneout, plotlssvm


% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

if isempty(varargin)
    kernel = 'RBF_kernel';
else
    kernel = varargin{1};
end

if type(1)=='f'
    perffun = 'mse';
elseif type(1)=='c'
    perffun = 'misclass';
else
    error('Type not supported. Choose ''f'' or ''c''')
end

n = size(x,1);
if n <= 300
    optfun = 'leaveoneoutlssvm';
    optargs = {perffun};
else
    optfun = 'crossvalidatelssvm';
    optargs = {10,perffun};
end

model = initlssvm(x,y,type,[],[],kernel);
model = tunelssvm(model,'simplex',optfun,optargs);
model = trainlssvm(model);

if size(x,2) <= 2
    plotlssvm(model);
end

% first output
yp = simlssvm(model,x);

% second output
alpha = model.alpha;

% third output
b = model.b;

% fourth and fifth output
gam = model.gam; sig2 = model.kernel_pars;


