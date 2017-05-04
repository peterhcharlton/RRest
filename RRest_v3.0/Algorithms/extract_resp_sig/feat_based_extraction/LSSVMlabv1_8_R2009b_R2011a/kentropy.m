function [en,U,lam] = kentropy(X, A1, A2, A3, A4)
% Quadratic Renyi Entropy for a kernel based estimator
% 
% Given the eigenvectors and the eigenvalues of the kernel matrix, the entropy is computed by
% 
% >> H = kentropy(X, U, lam)
% 
% The eigenvalue decomposition can also be computed (or
% approximated) implicitly:
% 
% >> H = kentropy(X, kernel, sig2)
% 
%
% Full syntax
% 
% >> H = kentropy(X, kernel, kernel_par)
% >> H = kentropy(X, kernel, kernel_par, type)
% >> H = kentropy(X, kernel, kernel_par, type, nb)
% 
%       Outputs    
%         H          : Quadratic Renyi entropy of the kernel matrix
%       Inputs    
%         X          : N x d matrix with the training data
%         kernel     : Kernel type (e.g. 'RBF_kernel')
%         kernel_par : Kernel parameter (bandwidth in the case of the 'RBF_kernel')
%         type(*)    : 'eig'(*), 'eigs', 'eign'
%         nb(*)      : Number of eigenvalues/eigenvectors used in the eigenvalue decomposition approximation
% 
%
% >> H = kentropy(X, U, lam)
% 
%       Outputs    
%         H   : Quadratic Renyi entropy of the kernel matrix
%       Inputs    
%         X   : N x d matrix with the training data
%         U   : N x nb matrix with principal eigenvectors
%         lam : nb x 1 vector with eigenvalues of principal components
% 
% See also:
%   kernel_matrix, RBF_kernel, demo_fixedsize

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

n= size(X,1);

if isstr(A1), % kernel_matrix

  kernel = A1;
  kernel_par = A2;
  eval('etype = A3;','etype =''eig'';');
  if ~(strcmp(etype, 'eig') |strcmp(etype, 'eigs') |strcmp(etype,'eign')),
    error('type has to be ''eig'', ''eigs'' or ''eign''...');
  end
  eval('nb = A4;',' ');
  
  if strcmp(etype,'eign'),
    eval('[U,lam] = eign(X,kernel,kernel_par,nb);','[U,lam] = eign(X,kernel,kernel_par);');
  else
    omega = kernel_matrix(X, kernel, kernel_par);
    eval('[U,lam] = feval(etype,omega,nb);','[U,lam] = feval(etype,omega);');
    if size(lam,1)==size(lam,2), lam = diag(lam); end
    %onen = ones(n,1)./n; en = -log(onen'*omega*onen);
  end
  

else
  U = A1;
  lam = A2;
end  
en = -log((sum(U,1)/n).^2 * lam);
  