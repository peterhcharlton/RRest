function cost = gcrossvalidate(model,estfct)

% Estimate the model performance of a model with genralized crossvalidation
% for regression with the LS-SVM
%
% >> cost = gcrossvalidate({Xtrain,Ytrain,type,gam,sig2})
% >> cost = gcrossvalidate( model)
% 
% Instead of dividing the data into $L$ disjoint sets, one takes the
% complete data and the effective degrees of freedom (effective number of parameters) 
% into account. The assumption is made that the input data are distributed 
% independent and identically over the input space. 
% 
% >> cost = gcrossvalidate(model)
% 
% Some commonly used criteria are:
% 
% >> cost = gcrossvalidate(model, 'misclass')
% >> cost = gcrossvalidate(model, 'mse')
% >> cost = gcrossvalidate(model, 'mae')
% 
% Full syntax
% 
%     1. Using LS-SVMlab with the functional interface:
% 
% >> cost = gcrossvalidate({X,Y,type,gam,sig2,kernel,preprocess}, estfct)
% 
%       Outputs    
%         cost          : Cost estimation of the generalized cross-validation
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
%
%     2. Using the object oriented interface:
% 
% >> cost = gcrossvalidate(model, estfct)
% 
%       Outputs    
%         cost          : Cost estimation of the generalized cross-validation
%
%       Inputs    
%         model         : Object oriented representation of the LS-SVM model
%         estfct(*)     : Function estimating the cost based on the residuals (by default mse)
% 
% 
%
% See also:
% leaveoneout, crossvalidate, trainlssvm, simlssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


% LS-SVMlab
eval('model = initlssvm(model{:});',' ');
eval('estfct;','estfct=''mse'';');

%
% initialisation and defaults
%
nb_data = size(model.ytrain,1);

% Y is raw data, non preprocessed
py = model.ytrain;
[~,Y] = postlssvm(model,[],py);

% Calculate kernel matrix and trace smoother
K = kernel_matrix(model.xtrain,model.kernel_type,model.kernel_pars);
tr = trace(smoother(model,K));

% Solve the linear system
S = ones(model.nb_data,1);
sol = linsolve([0 S';S K+eye(model.nb_data)./model.gam],[0;py],struct('SYM',true));

% Simulate
yh = K*sol(2:end) + ones(model.nb_data,1)*sol(1);
[~,yh] = postlssvm(model,[],yh);

% Generalized cross-validation
if ~(model.type(1)=='c')
    cost = feval(estfct,yh-Y);
else
    cost = feval(estfct,Y,sign(yh));
end
cost = cost/((1-tr/size(model.ytrain,1))^2);


function S = smoother(model,K)
% Smoother Matrix of the LS-SVM
% f = S*y
Z = pinv(K+eye(model.nb_data)./model.gam);
c = sum(sum(Z));
J = (ones(model.nb_data)./c);
S = K*(Z-Z*J*Z) + J*Z;



