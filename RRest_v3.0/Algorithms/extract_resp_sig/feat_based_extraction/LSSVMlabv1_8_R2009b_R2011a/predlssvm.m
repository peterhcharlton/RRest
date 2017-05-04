function pi = predlssvm(model,Xt,alpha,conftype)

% Construction of bias corrected 100(1-\alpha)% pointwise or
% simultaneous prediction intervals
%
% >> pi = predlssvm({X,Y,type,gam,kernel_par,kernel,preprocess}, Xt, alpha,conftype)
% >> pi = predlssvm(model, Xt, alpha, conftype)
%
% This function calculates bias corrected 100(1-\alpha)% pointwise or
% simultaneous prediction intervals. The procedure support homoscedastic
% data sets as well heteroscedastic data sets. The construction of the
% prediction intervals are based on the central limit theorem for linear
% smoothers combined with bias correction and variance estimation.
%
% 1. Using the functional interface:
%
%
% >> pi = predlssvm({X,Y,type,gam,kernel_par,kernel,preprocess})
% >> pi = predlssvm({X,Y,type,gam,kernel_par,kernel,preprocess}, alpha)
% >> pi = predlssvm({X,Y,type,gam,kernel_par,kernel,preprocess}, alpha, conftype)
%
%
%      Outputs
%        pi            : N x 2 matrix containing the lower and upper prediction intervals
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
%        conftype(*)   : Type of prediction interval 'pointwise' or 'simultaneous' (by default 'simultaneous')
%
% 2. Using the object oriented interface:
%
%
% >> pi = predlssvm(model)
% >> pi = predlssvm(model, alpha)
% >> pi = predlssvm(model, alpha, conftype)
%
%
%      Outputs
%        pi          : N x 2 matrix containing the lower and upper confidence intervals
%
%      Inputs
%        model       : Object oriented representation of the LS-SVM model
%        alpha       : Significance level (by default 5%)
%        conftype    : Type of prediction interval 'pointwise' or 'simultaneous' (by default 'simultaneous')
%
%
%  See also:
%    trainlssvm, simlssvm, cilssvm


% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

if iscell(model)
    model = initlssvm(model{:});
end

if nargin == 1
    error('Please specify a test set');
elseif nargin <= 2
    alpha = 0.05;
    conftype = 'simul';
    
elseif nargin <= 3
    conftype = 'simul';
end


if model.preprocess(1)=='p'
    [x,y] = postlssvm(model,model.xtrain,model.ytrain);
else
    x = model.xtrain; y = model.ytrain;
end

% train model
model = trainlssvm(model);

% smoother on test
S = smootherlssvm(model,Xt);

% bias: double smoothing with fourt order kernel RBF4
modelb = initlssvm(x,y,'f',[],[],'RBF4_kernel','o');
modelb = tunelssvm(modelb,'simplex','crossvalidatelssvm',{10,'mse'});
modelb = trainlssvm(modelb);

% output on test Yphat
Yphat = simlssvm(modelb,Xt);

% output of the smoother based on fourth order kernel on training data Yhat
Yhat = simlssvm(modelb,x);

biascorr = S*Yhat-Yphat;

%1) estimate variance nonparametrically
[sigma2,modele] = varest(model);
% smoother on training data
L = smootherlssvm(model);
sigma2t = varfun(modele,L,Xt);

%2) calculate var-cov matrix
S = S*diag(sigma2)*S';

%2b) find standardized absolute maxbias 
delta = max(abs(biascorr./sqrt(diag(S))));

%3) pointwise or simultaneous
if conftype(1)=='s'
    if model.preprocess(1)=='p'
        model.xtrain = prelssvm(model,x); 
    end
    z = tbform(model,alpha) + delta;
elseif conftype(1)=='p'
    z = norminv(alpha/2);
    Yphat = Yphat - biascorr;
else
    error('Wrong type of confidence interval. Please choose ''pointwise'' or ''simultaneous''');
end

V = sqrt(sigma2t+diag(S));
pi = [Yphat+z*V Yphat-z*V];

function V = varfun(modele,L,xt)
% estimate variance op new data, given the smooth of squared residuals
% "modele" and smoother matrix of the original smooth L
V = max(simlssvm(modele,xt),0);
% make estimate of V unbiased in homoscedastic case if regression
% estimate is unbiased
S = smootherlssvm(modele,xt);
V = V./(ones(size(xt,1),1)+S*diag(L*L'-2*L));

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

% make estimate of Var unbiased in homoscedastic case if regression
% estimate is unbiased
L = smootherlssvm(model);
S = smootherlssvm(modele);

var = var./(ones(size(x,1),1)+S*diag(L*L'-L-L'));