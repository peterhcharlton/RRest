function [sig_e, bay,model] = bay_errorbar(model,Xt, type, nb, bay)
% Compute the error bars for a one dimensional regression problem
% 
% >> sig_e = bay_errorbar({X,Y,'function',gam,sig2}, Xt)
% >> sig_e = bay_errorbar(model, Xt)
% 
% The computation takes into account the estimated noise variance
% and the uncertainty of the model parameters, estimated by
% Bayesian inference. sig_e is the estimated standard deviation of
% the error bars of the points Xt. A plot is obtained by replacing
% Xt by the string 'figure'.
% 
%
% Full syntax
% 
%     1. Using the functional interface:
% 
% >> sig_e = bay_errorbar({X,Y,'function',gam,sig2,kernel,preprocess}, Xt)
% >> sig_e = bay_errorbar({X,Y,'function',gam,sig2,kernel,preprocess}, Xt, type)
% >> sig_e = bay_errorbar({X,Y,'function',gam,sig2,kernel,preprocess}, Xt, type, nb)
% >> sig_e = bay_errorbar({X,Y,'function',gam,sig2,kernel,preprocess}, 'figure')
% >> sig_e = bay_errorbar({X,Y,'function',gam,sig2,kernel,preprocess}, 'figure', type)
% >> sig_e = bay_errorbar({X,Y,'function',gam,sig2,kernel,preprocess}, 'figure', type, nb)
% 
%       Outputs    
%         sig_e         : Nt x 1 vector with the [$ \sigma^2$] errorbands of the test data
%       Inputs    
%         X             : N x d matrix with the inputs of the training data
%         Y             : N x 1 vector with the inputs of the training data
%         type          : 'function estimation' ('f')
%         gam           : Regularization parameter
%         sig2          : Kernel parameter
%         kernel(*)     : Kernel type (by default 'RBF_kernel')
%         preprocess(*) : 'preprocess'(*) or 'original'
%         Xt            : Nt x d matrix with the inputs of the test data
%         type(*)       : 'svd'(*), 'eig', 'eigs' or 'eign'
%         nb(*)         : Number of eigenvalues/eigenvectors used in the eigenvalue decomposition approximation
%
%     2. Using the object oriented interface:
% 
% >> [sig_e, bay, model] = bay_errorbar(model, Xt)
% >> [sig_e, bay, model] = bay_errorbar(model, Xt,       type)
% >> [sig_e, bay, model] = bay_errorbar(model, Xt,       type, nb)
% >> [sig_e, bay, model] = bay_errorbar(model, 'figure')
% >> [sig_e, bay, model] = bay_errorbar(model, 'figure', type)
% >> [sig_e, bay, model] = bay_errorbar(model, 'figure', type, nb)
% 
%       Outputs    
%         sig_e     : Nt x 1 vector with the [$ \sigma^2$] errorbands of the test data
%         model(*)  : Object oriented representation of the LS-SVM model
%         bay(*)    : Object oriented representation of the results of the Bayesian inference
%       Inputs    
%         model     : Object oriented representation of the LS-SVM model
%         Xt        : Nt x d matrix with the inputs of the test data
%         type(*)   : 'svd'(*), 'eig', 'eigs' or 'eign'
%         nb(*)     : Number of eigenvalues/eigenvectors used in the eigenvalue decomposition approximation
% 
% See also:
%   bay_lssvm, bay_optimize, bay_modoutClass, plotlssvm

% Copyright (c) 2002,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.ac.be/sista/lssvmlab

if iscell(model), model = initlssvm(model{:}); end


if model.type(1)~='f',
  error(['confidence bounds only for function estimation. For' ...
	 ' classification, use ''bay_modoutClass(...)'' instead;']);
end

eval('type;','type=''svd'';');
eval('nb;','nb=model.nb_data;');
if ~(strcmpi(type,'svd') | strcmpi(type,'eig') | strcmpi(type,'eigs') | strcmpi(type,'eign')),
  error('Eigenvalue decomposition via ''svd'', ''eig'', ''eigs'' or ''eign''...');
end

if strcmpi(type,'eign')
  warning('The resulting errorbars are most probably not very usefull...');  
end

if ~isstr(Xt),

  eval('[sig_e, bay] = bay_confb(model,Xt,type,nb,bay);',...
       '[sig_e, bay] = bay_confb(model,Xt,type,nb);');
  
else

  grid = 50;
  [X,Y] = postlssvm(model,model.xtrain,model.ytrain);
  eval('[sig_e, bay] = bay_confb(model,X,type,nb,bay);',...
       '[sig_e, bay] = bay_confb(model,X,type,nb);');

  % plot the curve including confidence bound
  sige = sqrt(sig_e);
  Yt = simlssvm(model,X);
  
  figure;
  hold on;
  title(['LS-SVM_{\gamma=' num2str(model.gam(1)) ', \sigma^2=' num2str(model.kernel_pars(1)) ...
	 '}^{' model.kernel_type(1:3) '} and its 95% (2\sigma) error bands']);

  if model.x_dim==1,
    xlabel('X');
    ylabel('Y');
    [~,si] = sort(X);
    plot(X(si),Yt(si),'k'); hold on; 
    plot(X(si),Yt(si)+2.*sige(si),'-.r');    
    plot(X(si),Yt(si)-2.*sige(si),':r');
    plot(X(si),Y(si),'k*'); hold off;
  else
    xlabel('time');
    ylabel('Y');
    plot(Yt,'k'); hold on; 
    plot(Yt+2.*sige,'-.r'); 
    plot(Yt-2.*sige,':r');
    plot(Y,'k*'); hold off;
  end
 
end



function [sig_e, bay] = bay_confb(model,X,type,nb,bay)
% see formula's thesis TvG blz 126


nD = size(X,1);
%tol = .0001;

%
% calculate the eigenvalues
%
eval('bay;','[c1,c2,c3,bay] = bay_lssvm(model,1,type,nb);');
omega = kernel_matrix(model.xtrain(model.selector,1:model.x_dim), ...
		      model.kernel_type, model.kernel_pars);
oo = ones(1,model.nb_data)*omega;

% kernel values of  X
theta = kernel_matrix(model.xtrain(model.selector, 1:model.x_dim), ...
		      model.kernel_type, model.kernel_pars, X);
for i=1:nD,
  kxx(i,1) = feval(model.kernel_type, X(i,:),X(i,:), model.kernel_pars);
end



Zc = eye(model.nb_data) - ones(model.nb_data)./model.nb_data;


Hd = (Zc*bay.Rscores);
Hd = Hd*diag(1./bay.mu - (bay.mu+ bay.zeta*bay.eigvals).^-1)*Hd';


% forall x
for i=1:nD,
  term1(i,1) = bay.zeta^-1 + kxx(i)/bay.mu - theta(:,i)'*Hd*theta(:,i);
  term2(i,1) = 2/model.nb_data*sum(theta(:,i)'*Hd*omega) - 2/bay.mu/model.nb_data* sum(theta(:,i));
end


% once
term3 = 1/(bay.zeta*model.nb_data) ...
	+ 1/(bay.mu*model.nb_data^2)* sum(oo)  ...
	-1/(model.nb_data^2)* oo*Hd*oo';

sig_e = term1+term2+term3;


