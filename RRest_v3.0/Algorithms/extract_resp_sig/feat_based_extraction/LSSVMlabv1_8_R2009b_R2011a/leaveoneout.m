function cost = leaveoneout(model, estfct,combinefct)

% Estimate the performance of a trained model with leave-one-out crossvalidation
%
% CAUTION!! Use this function only to obtain the value of the leave-one-out score 
% function given the tuning parameters. Do not use this function together with 
% 'tunelssvm', but use 'leaveoneoutlssvm' instead. The latter is a faster 
% implementation which uses previously computed results.
%
% >> leaveoneout({X,Y,type,gam,sig2})
% >> leaveoneout(model)
%
% In each iteration, one leaves one point, and fits a model on the
% other data points. The performance of the model is estimated
% based on the point left out. This procedure is repeated for each
% data point. Finally, all the different estimates of the
% performance are combined (default by computing the mean). The
% assumption is made that the input data is distributed independent
% and identically over the input space.
%
%
% Full syntax
%
%     1. Using the functional interface for the LS-SVMs:
%
% >> cost = leaveoneout({X,Y,type,gam,sig2,kernel,preprocess})
% >> cost = leaveoneout({X,Y,type,gam,sig2,kernel,preprocess}, estfct)
% >> cost = leaveoneout({X,Y,type,gam,sig2,kernel,preprocess}, estfct, combinefct)
%
%       Outputs
%         cost          : Cost estimated by leave-one-out crossvalidation
%
%       Inputs
%         X             : Training input data used for defining the LS-SVM and the preprocessing
%         Y             : Training output data used for defining the LS-SVM and the preprocessing
%         type          : 'function estimation' ('f') or 'classifier' ('c')
%         gam           : Regularization parameter
%         sig2          : Kernel parameter (bandwidth in the case of the 'RBF_kernel')
%         kernel(*)     : Kernel type (by default 'RBF_kernel')
%         preprocess(*) : 'preprocess'(*) or 'original'
%         estfct(*)     : Function estimating the cost based on the residuals (by default mse)
%         combinefct(*) : Function combining the estimated costs on the different folds (by default mean)
%
%
%     2. Using the object oriented interface for the LS-SVMs:
%
% >> cost = leaveoneout(model)
% >> cost = leaveoneout(model, estfct)
% >> cost = leaveoneout(model, estfct, combinefct)
%
%       Outputs
%         cost          : Cost estimated by leave-one-out crossvalidation
%
%       Inputs
%         model         : Object oriented representation of the model
%         estfct(*)     : Function estimating the cost based on the residuals (by default mse)
%         combinefct(*) : Function combining the estimated costs on the different folds (by default mean)
%
%
%
% See also:
%    crossvalidate, trainlssvm, simlssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

% LS-SVMlab
eval('model = initlssvm(model{:});',' ');
eval('estfct;','estfct=''mse'';');
eval('combinefct;','combinefct=''mean'';');

%
%initialize: no incremental  memory allocation
%
p = randperm(model.nb_data); 
px = model.xtrain(p,:);
py = model.ytrain(p,:);
[~,Y] = postlssvm(model,[],py); % Y is raw data, non preprocessed

% kernel matrix computation
K = kernel_matrix(px,model.kernel_type,model.kernel_pars);

Ka = pinv([K+eye(model.nb_data)./model.gam ones(model.nb_data,1);ones(1,model.nb_data) 0]);
sol = Ka*[py;0]; model.alpha = sol(1:end-1); model.b = sol(end);
yh = py - model.alpha./diag(Ka(1:model.nb_data,1:model.nb_data));

[~,yh] = postlssvm(model,[],yh);
if ~(model.type(1)=='c')
    cost = feval(estfct,yh-Y);
else
    cost = feval(estfct,Y,sign(yh));
end


