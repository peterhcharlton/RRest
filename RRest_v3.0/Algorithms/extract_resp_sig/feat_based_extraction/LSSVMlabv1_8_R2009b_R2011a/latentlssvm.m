function [zt,model] = latentlssvm(varargin)
% Calculate the latent variables of the LS-SVM classifier at the given test data
% 
% >> Zt = latentlssvm({X,Y,'classifier',gam,sig2,kernel}, {alpha,b}, Xt)
% >> Zt = latentlssvm({X,Y,'classifier',gam,sig2,kernel}, Xt)
% >> [Zt, model] = latentlssvm(model, Xt)
% 
% The latent variables of a binary classifier are the continuous
% simulated values of the test data which are used to make the
% final classifications. The classification of a testpoint depends
% on whether the latent value exceeds the model's threshold (b). If
% appropriate, the model is trained by the standard procedure (trainlssvm) first.
% 
% As an application example: crossvalidation can be based on the latent variables:
% 
% >> cost = crossvalidate(model, X, Y, 10, 'mse', 'mean', 'original', 'trainlssvm', 'latentlssvm')
% 
%
% Full syntax
% 
%     1. Using the functional interface:
% 
% >> Zt = latentlssvm({X,Y,type,gam,sig2,kernel,preprocess}, Xt)
% 
%       Outputs    
%         Zt            : Nt x m matrix with predicted latent simulated outputs
%       Inputs    
%         X             : N x d matrix with the inputs of the training data
%         Y             : N x 1 vector with the outputs of the training data
%         type          : 'classifier' ('c')
%         gam           : Regularization parameter
%         sig2          : Kernel parameter (bandwidth in the case of the 'RBF_kernel')
%         kernel(*)     : Kernel type (by default 'RBF_kernel')
%         preprocess(*) : 'preprocess'(*) or 'original'
%         Xt            : Nt x d matrix with the inputs of the test data
% 
%
%     2. Using the object oriented interface:
% 
% >> [Zt, model] = latentlssvm(model, Xt)
% 
%       Outputs    
%         Zt       : Nt x m matrix with continuous latent simulated outputs
%         model(*) : Trained object oriented representation of the LS-SVM model
%       Inputs    
%         model    : Object oriented representation of the LS-SVM model
%         Xt       : Nt x d matrix with the inputs of the test data
% 
% See also:
%   trainlssvm, simlssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

model = varargin{1};
if iscell(model),
  model = initlssvm(model{:});
end

if model.type(1)~='c',
  error('Only usefull for classification tasks...');
end
[~, zt, model] = simlssvm(varargin{:});