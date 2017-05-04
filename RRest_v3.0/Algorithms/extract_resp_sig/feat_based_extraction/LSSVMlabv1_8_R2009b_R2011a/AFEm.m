function [features,eigvec,eigvals] = AFEm(Xs,kernel, kernel_pars,X,type,nb,eigvec,eigvals)
% Automatic Feature Extraction by Nystrom method
%
%
% >> features = AFE(X, kernel, sig2, Xt)
%
% Description
% Using the Nystr�m approximation method, the mapping of data to
% the feature space can be evaluated explicitly. This gives the
% features that one can use for a linear regression or
% classification. The decomposition of the mapping to the feature
% space relies on the eigenvalue decomposition of the kernel
% matrix. The Matlab ('eigs') or Nystr�m's ('eign') approximation
% using the nb most important eigenvectors/eigenvalues can be
% used. The eigenvalue decomposition is not re-calculated if it is
% passed as an extra argument. This routine internally calls a cmex file.
%
% Full syntax
% 
% >> [features, U, lam] = AFE(X, kernel, sig2, Xt) 
% >> [features, U, lam] = AFE(X, kernel, sig2, Xt, type) 
% >> [features, U, lam] = AFE(X, kernel, sig2, Xt, type, nb) 
% >> features          = AFE(X, kernel, sig2, Xt, [],[], U, lam)
% 
% Outputs    
%   features : Nt x nb matrix with extracted features
%   U(*)     : N x nb matrix with eigenvectors
%   lam(*)   : nb x 1 vector with eigenvalues
% Inputs    
%   X      : N x d matrix with input data
%   kernel : Name of the used kernel (e.g. 'RBF_kernel')
%   sig2   : parameter of the used kernel
%   Xt     : Data from which the features are extracted
%   type(*): 'eig'(*), 'eigs' or 'eign'
%   nb(*)  : Number of eigenvalues/eigenvectors used in the eigenvalue decomposition approximation
%   U(*)   : N x nb matrix with eigenvectors
%   lam(*) : nb x 1 vector with eigenvalues
% 
% See also:
%   kernel_matrix, RBF_kernel, demo_fixedsize

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

N = size(X,1);
Nc = size(Xs,1);

eval('type;','type=''eig'';');
if ~(strcmp(type,'eig') || strcmp(type,'eigs') || strcmp(type,'eign') )
  error('Type needs to be ''eig'', ''eigs'' or ''eign''...');
end
  

% eigenvalue decomposition to do..
if nargin<=7,
  omega = kernel_matrix(Xs, kernel, kernel_pars);
  if strcmp(type,'eig'),
    [eigvec,eigvals] = eig(omega+2*eye(size(omega,1))); % + jitter factor
    eigvals = diag(eigvals); 
    clear omega
  elseif strcmp(type,'eigs'),
    eval('nb;','nb=min(size(omega,1),10);'); options.disp = 0;
    [eigvec,eigvals] = eigs(omega+2*eye(size(omega,1)),nb,'lm',options); clear omega % + jitter factor
  elseif strcmp(type,'eign'),
    eval('nb;','nb=min(size(omega,1),10);'); 
    [eigvec,eigvals] = eign(omega+2*eye(size(omega,1)),nb); clear omega % + jitter factor
  end
  eigvals = (eigvals-2)/Nc;

  peff = eigvals>eps;
  eigvals = eigvals(peff);
  eigvec = eigvec(:,peff); clear peff
  
end 

if strcmp(kernel,'RBF_kernel')
    omegaN = sum(X.^2,2)*ones(1,Nc);
    omegaN = omegaN + ones(N,1)*sum(Xs.^2,2)';
    omegaN = omegaN -2*X*Xs'; clear X Xs
    omegaN = exp(-omegaN./kernel_pars);
    
elseif strcmp(kernel,'lin_kernel')
    omegaN = X*Xs'; clear X Xs
    
elseif strcmp(kernel,'poly_kernel')
    omegaN =  X*Xs'; clear X Xs
    omegaN = (omegaN + kernel_pars(1)).^kernel_pars(2);
else
    disp('kernel_type unkown')
    return;
end

%A=sqrt(Nc) ./ sqrt(eigvals); clear eigvals
features = omegaN*eigvec; clear omegaN
features = repmat((1 ./ sqrt(eigvals))',N,1).*features;


