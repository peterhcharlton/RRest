function [model, ss] = bay_initlssvm(model)
% Initialize the hyperparameters [$ \gamma$]  and [$ \sigma^2$]  before optimization with bay_optimize
% 
% >> [gam, sig2] = bay_initlssvm({X,Y,type,[],[]})
% >> model       = bay_initlssvm(model)
% 
% 
% A starting value for sig2 is only given if the model has kernel type 'RBF_kernel'.
% 
%
% Full syntax
% 
%     1. Using the functional interface:
% 
% >> [gam, sig2] = bay_initlssvm({X,Y,type,[],[],kernel})
% 
%       Outputs    
%         gam  : Proposed initial regularization parameter
%         sig2 : Proposed initial 'RBF_kernel' parameter
%       Inputs    
%         X    : N x d matrix with the inputs of the training data
%         Y    : N x 1 vector with the outputs of the training data
%         type : 'function estimation' ('f') or 'classifier' ('c')
%         kernel(*) : Kernel type (by default 'RBF_kernel')
%
%     2. Using the object oriented interface:
% 
% >> model = bay_initlssvm(model)
% 
%       Outputs    
%         model : Object oriented representation of the LS-SVM model with initial hyperparameters
%       Inputs    
%         model : Object oriented representation of the LS-SVM model
% 
% See also:
%    bay_lssvm, bay_optimize

% -disclaimer


if iscell(model),
  iscell_model = 1;
  model = initlssvm(model{:});
else
  iscell_model = 0;
end

% start sig2
% sig2 as the std of the fitting Gaussian
if strcmp(model.kernel_type,'RBF_kernel'),
  model.kernel_pars = sum(range(model.xtrain))./1.96.*ones(1,model.y_dim);
else
  warning('Only usefull for ''RBF_kernel''...');
  if iscell_model,
    ss = model.kernel_pars;
    model = model.gam;
  end
  return
end

% set starting value
if numel(model.gam)~=1, model.gam = 1; end


% start gamma
for i=1:10,
  gam(i,1) = exp(i-5);
  model.gam = gam(i,1);
  [~,~,~,bay] = bay_lssvm(model,2,'svd');
  gam(i,2) = bay.Geff;
end
[~,index] = sort(abs(gam(:,1)-gam(:,2).*2));
model.gam = gam(index(1),1);

if iscell_model,
  ss = model.kernel_pars;
  model = model.gam;
end
