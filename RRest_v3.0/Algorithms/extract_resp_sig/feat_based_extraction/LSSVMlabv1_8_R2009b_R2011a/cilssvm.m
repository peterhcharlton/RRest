function ci = cilssvm(model,alpha,conftype)
%
% Construction of bias corrected 100(1-\alpha)% pointwise or
% simultaneous confidence intervals
%
% >> ci = cilssvm({X,Y,type,gam,kernel_par,kernel,preprocess},alpha,conftype)
% >> ci = cilssvm(model,alpha,conftype)
%
% This function calculates bias corrected 100(1-\alpha)% pointwise or 
% simultaneous confidence intervals. The procedure support homoscedastic 
% data sets as well heteroscedastic data sets. The construction of the  
% confidence intervals are based on the central limit theorem for linear  
% smoothers combined with bias correction and variance estimation. 
%
% 1. Using the functional interface:
%
%
% >> ci = cilssvm({X,Y,type,gam,kernel_par,kernel,preprocess})
% >> ci = cilssvm({X,Y,type,gam,kernel_par,kernel,preprocess}, alpha)
% >> ci = cilssvm({X,Y,type,gam,kernel_par,kernel,preprocess}, alpha, conftype)
%
%
%      Outputs
%        ci            : N x 2 matrix containing the lower and upper confidence intervals
%  
%      Inputs
%        X             : N x d matrix with the inputs of the training data
%        Y             : N x 1 vector with the outputs of the training data
%        type          : 'function estimation' ('f') or 'classifier' ('c')
%        gam           : Regularization parameter
%        sig2          : Kernel parameter(s) (bandwidth in the case of the 'RBF_kernel')
%        kernel(*)     : Kernel type (by default 'RBF_kernel')
%        preprocess(*) : 'preprocess'(*) or 'original'
%        alpha(*)      : Significance level (by default 5%)
%        conftype(*)   : Type of confidence interval 'pointwise' or 'simultaneous' (by default 'simultaneous')
%
% 2. Using the object oriented interface:
%
%
% >> ci = cilssvm(model)
% >> ci = cilssvm(model, alpha)
% >> ci = cilssvm(model, alpha, conftype)
%
%
%      Outputs
%        ci          : N x 2 matrix containing the lower and upper confidence intervals
%  
%      Inputs
%        model       : Object oriented representation of the LS-SVM model
%        alpha       : Significance level (by default 5%)
%        conftype    : Type of confidence interval 'pointwise' or 'simultaneous' (by default 'simultaneous')
%
%
%  See also:
%    trainlssvm, simlssvm, predlssvm


% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

if iscell(model)
    model = initlssvm(model{:});
end

if nargin <= 1
    alpha = 0.05;
    conftype = 'simul';
    
elseif nargin <= 2
    conftype = 'simul';
end

if model.preprocess(1)=='p'
    error('Please use original data to compute confidence intervals...')
end

x = model.xtrain; y = model.ytrain;

% train model
if isempty(model.gam) && isempty(model.kernel.pars)
    error('Please tune model first with ''tunelssvm'' to obtain tuning parameters');
end
model = trainlssvm(model);

s = smootherlssvm(model);
Yhat = simlssvm(model,x);

% bias: double smoothing with fourt order kernel RBF4
modelb = initlssvm(x,y,'f',[],[],'RBF4_kernel','o');
modelb = tunelssvm(modelb,'simplex','crossvalidatelssvm',{10,'mse'});
modelb = trainlssvm(modelb);

biascorr = (s-eye(size(x,1)))*simlssvm(modelb,x);

% construct approximate 100(1-alpha)% confidence interval
%1) estimate variance nonparametrically
sigma2 = varest(model);

%2) calculate var-cov matrix
s = s*diag(sigma2)*s';

%2b) find standardized absolute maxbias 
delta = max(abs(biascorr./sqrt(diag(s))));

%3) pointwise or simultaneous?
if conftype(1)=='s'
    z = tbform(model,alpha) + delta;
elseif conftype(1)=='p'
    z = norminv(alpha/2);
    Yhat = Yhat - biascorr;
else
    error('Wrong type of confidence interval. Please choose ''pointwise'' or ''simultaneous''');
end
    
ci = [Yhat+z*sqrt(diag(s)) Yhat-z*sqrt(diag(s))];

function [var,modele] = varest(model)

% if preprocessed data, construct original data
if model.preprocess(1)=='p'
    [x,y] = postlssvm(model,model.xtrain,model.ytrain);
else
    x = model.xtrain; y = model.ytrain;
end

model = trainlssvm(model);

Yh = simlssvm(model,x);

% Squared normalized residuals
e2 = (y-Yh).^2;

% Make variance model
if model.nb_data <= 200
    costfun = 'leaveoneoutlssvm'; costargs = {'mae'};
else
    costfun = 'crossvalidatelssvm'; costargs = {10,'mae'};
end
modele = initlssvm(x,e2,'f',[],[],'RBF_kernel');
modele = tunelssvm(modele,'simplex',costfun,costargs);
modele = trainlssvm(modele);

% variance model
var = max(simlssvm(modele,x),0);

% make estimate of var unbiased in homoscedastic case if regression
% estimate is unbiased
L = smootherlssvm(model);
S = smootherlssvm(modele);

var = var./(ones(size(x,1),1)+S*diag(L*L'-L-L'));