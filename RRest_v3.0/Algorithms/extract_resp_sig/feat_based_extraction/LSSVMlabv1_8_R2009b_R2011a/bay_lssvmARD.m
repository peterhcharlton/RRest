function [inputs,ordered,costs,sig2n,model] = bay_lssvmARD(model,type,btype,nb);
% Bayesian Automatic Relevance Determination of the inputs of an LS-SVM
% 
% 
% >> dimensions = bay_lssvmARD({X,Y,type,gam,sig2})
% >> dimensions = bay_lssvmARD(model)
% 
% For a given problem, one can determine the most relevant inputs
% for the LS-SVM within the Bayesian evidence framework. To do so,
% one assigns a different weighting parameter to each dimension in
% the kernel and optimizes this using the third level of
% inference. According to the used kernel, one can remove inputs
% corresponding the larger or smaller kernel parameters. This
% routine only works with the 'RBF_kernel' with a sig2 per
% input. In each step, the input with the largest optimal sig2 is
% removed (backward selection). For every step, the generalization
% performance is approximated by the cost associated with the third
% level of Bayesian inference.
% 
% The ARD is based on backward selection of the inputs based on the
% sig2s corresponding in each step with a minimal cost
% criterion. Minimizing this criterion can be done by 'continuous'
% or by 'discrete'. The former uses in each step continuous varying
% kernel parameter optimization, the latter decides which one to
% remove in each step by binary variables for each component (this
% can only applied for rather low dimensional inputs as the number
% of possible combinations grows exponentially with the number of
% inputs). If working with the 'RBF_kernel', the kernel parameter
% is rescaled appropriately after removing an input variable.
% 
% The computation of the Bayesian cost criterion can be based on
% the singular value decomposition 'svd' of the full kernel matrix
% or by an approximation of these eigenvalues and vectors by the
% 'eigs' or 'eign' approximation based on 'nb' data points.
% 
% Full syntax
% 
%     1. Using the functional interface:
% 
% >> [dimensions, ordered, costs, sig2s] =  bay_lssvmARD({X,Y,type,gam,sig2,kernel,preprocess})
% >> [dimensions, ordered, costs, sig2s] =  bay_lssvmARD({X,Y,type,gam,sig2,kernel,preprocess}, method)
% >> [dimensions, ordered, costs, sig2s] =  bay_lssvmARD({X,Y,type,gam,sig2,kernel,preprocess}, method, type)
% >> [dimensions, ordered, costs, sig2s] =  bay_lssvmARD({X,Y,type,gam,sig2,kernel,preprocess}, method, type, nb)
% 
%       Outputs    
%         dimensions : r x 1 vector of the relevant inputs
%         ordered(*) : d x 1 vector with inputs in decreasing order of relevance
%         costs(*)   : Costs associated with third level of inference in every selection step
%         sig2s(*)   : Optimal kernel parameters in each selection step
%       Inputs    
%         X          : N x d matrix with the inputs of the training data
%         Y          : N x 1 vector with the outputs of the training data
%         type       : 'function estimation' ('f') or 'classifier' ('c')
%         gam        : Regularization parameter
%         sig2       : Kernel parameter (bandwidth in the case of the 'RBF_kernel')
%         kernel(*)  : Kernel type (by default 'RBF_kernel')
%         preprocess(*) :'preprocess'(*) or 'original'
%         method(*)  : 'discrete'(*) or 'continuous'
%         type(*)    : 'svd'(*), 'eig', 'eigs', 'eign'
%         nb(*)      :Number of eigenvalues/eigenvectors used in the eigenvalue decomposition approximation
% 
%     2. Using the object oriented interface:
% 
% >> [dimensions, ordered, costs, sig2s, model] = bay_lssvmARD(model)
% >> [dimensions, ordered, costs, sig2s, model] = bay_lssvmARD(model, method)
% >> [dimensions, ordered, costs, sig2s, model] = bay_lssvmARD(model, method, type)
% >> [dimensions, ordered, costs, sig2s, model] = bay_lssvmARD(model, method, type, nb)
% 
%       Outputs    
%         dimensions : r x 1 vector of the relevant inputs
%         ordered(*) : d x 1 vector with inputs in decreasing order of relevance
%         costs(*)   : Costs associated with third level of inference in every selection step
%         sig2s(*)   : Optimal kernel parameters in each selection step
%         model(*)   : Object oriented representation of the LS-SVM model trained only on the relevant inputs
%       Inputs    
%         model      : Object oriented representation of the LS-SVM model
%         method(*)  : 'discrete'(*) or 'continuous'
%         type(*)    : 'svd'(*), 'eig', 'eigs', 'eign'
%         nb(*)      : Number of eigenvalues/eigenvectors used in the eigenvalue decomposition approximation
% 
% See also:
%   bay_lssvm, bay_optimize, bay_modoutClass, bay_errorbar


% Copyright (c) 2002,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.ac.be/sista/lssvmlab




warning off
eval('type;','type=''discrete'';');
eval('btype;','btype=''svd'';');

if ~(type(1)=='d' | type(1)=='c'),
  error('type needs to be ''continuous'' or ''discrete''...');
end

  
if ~(strcmpi(btype,'svd') | strcmpi(btype,'eig') | strcmpi(btype,'eigs') | strcmpi(btype,'eign')),
  error('Eigenvalue decomposition via ''svd'', ''eig'', ''eigs'' or ''eign''.');
end

% 
% initiate model
%
if ~isstruct(model), 
  model = initlssvm(model{:}); 
end

'OPTIMIZING GAMMA AND KERNEL PARAMETERS WITH BAYESIAN FRAMEWORK OVER ALL INPUTS...'

%model = changelssvm(model, 'kernel_type', 'RBF_kernel');
eval('[model,kernel_pars,bay] = bay_optimize(model,3,btype,nb);',...
     '[model,kernel_pars,bay] = bay_optimize(model,3,btype);');

costs(1) = bay.costL3;

%
% init parameters
%
eval('nb;','nb=inf;');
xdim = model.x_dim;
all = 1:xdim;
reject = zeros(xdim,1);




%
% continuous
%
if type(1)=='c',

  if length(model.kernel_pars)~=model.x_dim,
    model = changelssvm(model,'kernel_pars',model.kernel_pars(1)*ones(1,model.x_dim));
  end

  sig2n = zeros(xdim-1,xdim);
  [Xtrain,Ytrain] = postlssvm(model,model.xtrain(model.selector,:),model.ytrain(model.selector,:));
  for d=1:xdim-1,
    ['testing for ' num2str(xdim-d+1) ' inputs']
    [modelf,sig2n(d,1:(xdim-d+1)),bay] = ...
	bay_optimize({Xtrain(:,all), Ytrain,model.type,model.gam, model.kernel_pars(:,all),model.kernel_type, model.preprocess},...
		     3,btype,nb)
    costs(d+1,:) = bay.costL3;
    [m,reject(d)] = max(sig2n(d,:));
    all = setdiff(all,reject(d)); all=reshape(all,1,length(all));
    
    ['SELECTED INPUT(S) (''continuous''): [' num2str(all) ']']
  end
  reject(xdim) = all;
  
  
%
% discrete 
%
elseif type(1)=='d',   

  
  if length(model.kernel_pars)>1,
    error('only 1 kernel parameter supported for the moment, use ''fmin'' instead;');
  end
  [Xtrain,Ytrain] = postlssvm(model,model.xtrain(model.selector,:), ...
			      model.ytrain(model.selector,:));

  
  %
  % cost for all
  % 
  [c3,bay] = bay_lssvm({Xtrain, Ytrain,...
		    model.type,model.gam, model.kernel_pars,...
		    model.kernel_type, model.preprocess}, 3,btype,nb);    
  costs(1,:) = bay.costL3;
    
  
  %
  % iteration
  %
  for d=1:xdim-1,

    ['testing for ' num2str(xdim-d+1) ' inputs']
    
    % rescaling of kernel parameters
    if strcmp(model.kernel_type,'RBF_kernel'), 
      % RBF
      model.kernel_pars = (model.x_dim-d)./model.x_dim.*model.kernel_pars;
    else
      % else
      model = bay_optimize({Xtrain(:,all), Ytrain,...
		    model.type,model.gam, model.kernel_pars,model.kernel_type, model.preprocess},3,btype,nb);      
    end

    % which input to remove?
    minc3 = inf;
    for a = 1:length(all),
      [c3,bayf] = bay_lssvm({Xtrain(:,all([1:a-1 a+1:end])), Ytrain,...
		      model.type,model.gam, model.kernel_pars,...
		      model.kernel_type, model.preprocess}, 3,btype,nb);
      if c3<minc3, 
	minc3=c3;reject(d)=all(a);bay=bayf; 
      end
    end
    
    % remove input d...
    all = setdiff(all,reject(d));all=reshape(all,1,length(all));
    costs(d+1) = bay.costL3;
    %save ARD_ff
    
    ['SELECTED INPUT(S) (''discrete''): [' num2str(all) ']']

  end
  reject(xdim) = all;
  
end


%
% select best reduction (costL2 lowest)
%
[mcL2,best] = min(costs);
ordered = reject(end:-1:1);
inputs = ordered(1:xdim-(best-1));
eval('mkp = model.kernel_pars(:,inputs);','mkp = model.kernel_pars;');
model = initlssvm(Xtrain(:,inputs),Ytrain,model.type,model.gam, mkp, model.kernel_type, model.preprocess);

warning on
      