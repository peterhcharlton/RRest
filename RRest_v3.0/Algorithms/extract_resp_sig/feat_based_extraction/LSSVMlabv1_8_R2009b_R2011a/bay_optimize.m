function [model,A,B,C,D] = bay_optimize(model,level, type, nb, bay)
% Optimize the posterior probabilities of model (hyper-) parameters with respect to the different levels in Bayesian inference
% 
% One can optimize on the three different inference levels:
% 
%     - First level: In the first level one optimizes the support values alpha 's and the bias b.
%     - Second level: In the second level one optimizes the regularization parameter gam.
%     - Third level: In the third level one optimizes the kernel
%                    parameter. In the case of the common 'RBF_kernel' the kernel
%                    parameter is the bandwidth sig2. 
% This routine is only tested with Matlab version 6 using the corresponding optimization toolbox.
% 
% Full syntax
% 
%     1. Outputs on the first level:
% 
% >> [model, alpha, b] = bay_optimize({X,Y,type,gam,sig2,kernel,preprocess}, 1)
% >> [model, alpha, b] = bay_optimize(model, 1)
% 
%         model    : Object oriented representation of the LS-SVM model optimized on the first level
%         alpha(*) : Support values optimized on the first level of inference
%         b(*)     : Bias term optimized on the first level of inference
% 
%
%     2. Outputs on the second level:
% 
% >> [model,gam] = bay_optimize({X,Y,type,gam,sig2,kernel,preprocess}, 2)
% >> [model,gam] = bay_optimize(model, 2)
% 
%       model   : Object oriented representation of the LS-SVM model optimized on the second level of inference
%       gam(*)  : Regularization parameter optimized on the second level of inference
% 
%
%     3. Outputs on the third level:
% 
% >> [model, sig2] = bay_optimize({X,Y,type,gam,sig2,kernel,preprocess}, 3)
% 
%       model   : Object oriented representation of the LS-SVM model optimized on the third level of inference
%       sig2(*) : Kernel parameter optimized on the third level of inference
% 
%
%     4. Inputs using the functional interface
% 
% >> model = bay_optimize({X,Y,type,gam,sig2,kernel,preprocess}, level)
% >> model = bay_optimize({X,Y,type,gam,sig2,kernel,preprocess}, level, type)
% >> model = bay_optimize({X,Y,type,gam,sig2,kernel,preprocess}, level, type, nb)
% 
%         X             : N x d matrix with the inputs of the training data
%         Y             : N x 1 vector with the outputs of the training data
%         type          : 'function estimation' ('f') or 'classifier' ('c')
%         gam           : Regularization parameter
%         sig2          : Kernel parameter (bandwidth in the case of the 'RBF_kernel')
%         kernel(*)     : Kernel type (by default 'RBF_kernel')
%         preprocess(*) : 'preprocess'(*) or 'original'
%         level         : 1, 2, 3
%         type(*)       : 'eig', 'svd'(*), 'eigs', 'eign'
%         nb(*)         : Number of eigenvalues/eigenvectors used in the eigenvalue decomposition approximation
% 
%
%     5. Inputs using the object oriented interface
% 
% >> model = bay_optimize(model, level)
% >> model = bay_optimize(model, level, type)
% >> model = bay_optimize(model, level, type, nb)
% 
%         model   : Object oriented representation of the LS-SVM model
%         level   : 1, 2, 3
%         type(*) : 'eig', 'svd'(*), 'eigs', 'eign'
%         nb(*)   : Number of eigenvalues/eigenvectors used in the eigenvalue decomposition approximation
% 
% See also:
%   bay_lssvm, bay_lssvmARD, bay_modoutClass, bay_errorbar


% Copyright (c) 2002,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.ac.be/sista/lssvmlab

vers = version;
if vers(1)<'6', 
  error(['This routine is only supported currently under MATLAB 6' ...
	 ' and its corresponding optimization toolbox.']);
end

if iscell(model), model = initlssvm(model{:}); end


if ~(level==1 | level==2 | level==3),
  error('level must be 1, 2 or 3.');
end

eval('nb;','nb=model.nb_data;');



if level==1,

  eval('type;','type=''train'';');
  eval('[F1, F2, F3, C , model] = bay_lssvm(model,1,type,nb,bay);',...
       '[F1, F2, F3, C , model] = bay_lssvm(model,1,type,nb);');
  A = model.alpha
  B = model.b;
  
elseif level==2,

  eval('type;','type=''svd'';');
  eval('[model, A,B] = bay_optimize2(model,type,nb,bay); ',...
       '[model, A,B] = bay_optimize2(model,type,nb);')
  
elseif level==3,
  
  
  % check fminunc
  resp = which('fminunc');
  disp(' ');
  if isempty(resp),
    error(' ''fminunc'' not available');
  end

  eval('type;','type=''svd'';');
  
  %startvalues
  model = bay_optimize2(model,type,nb);
  % start value given in model: not fixed 'cause updating
  % of optimal parameters needs to be possible
  start_param = model.kernel_pars;

  opties=optimset('MaxFunEvals', 250, 'TolFun', .001, 'TolX', .001 );
  eval('A = fminunc(@costL3, start_param, opties, model, type, nb);');
 
  model = changelssvm(model,'kernel_pars',abs(A)); 
  [~,B,model] = bay_lssvm(model,3, type, nb);
end





function [model, A,B] = bay_optimize2(model,type,nb,bay)  

  % check fminunc
  resp = which('fminunc');
  disp(' ');
  if isempty(resp),
    error(' ''fminunc'' not available');
  end

  opties=optimset('TypicalX',model.kernel_pars,'MaxFunEvals', 2000,'GradObj','on','DerivativeCheck', 'off', 'TolFun', .0001, 'TolX', .0001 );
  if nargin<4,
    [c,dc,o, bay] = bay_lssvm(model,2,type,nb);
  end  
  eval('gam_opt = fminunc(@costL2, abs(model.gam), opties, model, type, nb,bay);');
  model = changelssvm(model,'gam',abs(gam_opt));
  [D1, D2, D3,B,model] = bay_lssvm(model,2,type, nb, bay);
  A = model.gam;

  

function [cost,Dcost] = costL2(lgam, model, type, nb, bay)
%
  model = changelssvm(model,'gam',abs(lgam+1000*eps));
  [cost, Dcost] = bay_lssvm(model,2,type, nb, bay);


function cost = costL3(sig2, model, type, nb)
%
  model = changelssvm(model,'kernel_pars',abs(sig2));
  cost = bay_lssvm(model,3, type, nb); 
  disp(['sig2 = ' num2str(model.kernel_pars) ' costL3 = ' num2str(cost) ';'])
