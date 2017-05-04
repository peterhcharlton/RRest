function [Xd,lam,U] = denoise_kpca(Xo,A1,A2,A3,A4,A5)
% Reconstruct the data mapped on the first principal components
% 
% >> Xd = denoise_kpca(X, kernel, kernel_par);
% 
% Denoising can be done by moving the point in inputspace so that
% its corresponding map to feature space is optimized. This means
% that the data point in feature space is as close as possible with
% its corresponding reconstructed points using the principal
% components. If the principal components are to be calculated on
% the same data X as the one one wants to denoise, use the command:
% 
% >> Xd         = denoise_kpca(X, kernel, kernel_par);
% >> [Xd,lam,U] = denoise_kpca(X, kernel, kernel_par, [], type, nb);
% 
% When one wants to denoise data 'Xt' other than the data used to obtain the principal components:
% 
% >> Xd          = denoise_kpca(X, kernel, kernel_par, Xt);
% >> [Xd, lam, U] = denoise_kpca(X, kernel, kernel_par, Xt, type, nb);
% 
%
% Full syntax
% 
% >> [Xd, lam, U] = denoise_kpca(X, kernel, kernel_par, Xt);
% >> [Xd, lam, U] = denoise_kpca(X, kernel, kernel_par, Xt, type);
% 
%       Outputs    
%         Xd         : N x d (Nt x d) matrix with denoised data X (Xt)
%         lam(*)     : nb x 1 vector with eigenvalues of principal components
%         U(*)       : N x nb (Nt x d) matrix with principal eigenvectors
%       Inputs    
%         X          : N x d matrix with data points used for finding the principal components
%         kernel     : Kernel type (e.g. 'RBF_kernel')
%         kernel_par : Kernel parameter (bandwidth in the case of the 'RBF_kernel')
%         Xt(*)      : Nt x d matrix with the points to denoise (if not specified, X is denoised instead)
%         type(*)    : 'eig'(*), 'svd', 'eigs', 'eign'
%         nb(*)      : Number of principal components used in approximation
% 
% >> Xd  = denoise_kpca(X, U, lam, kernel, kernel_par, Xt);
% 
%       Outputs    
%         Xd         : N x d (Nt x d) matrix with denoised data X (Xt)
%       Inputs    
%         X          : N x d matrix with data points used for finding the principal components
%         U          : N x nb (Nt x d) matrix with principal eigenvectors
%         lam        : nb x 1 vector with eigenvalues of principal components
%         kernel     : Kernel type (e.g. 'RBF_kernel')
%         kernel_par : Kernel parameter (bandwidth in the case of the 'RBF_kernel')
%         Xt(*)      : Nt x d matrix with the points to denoise (if not specified, X is denoised instead)
% 
% See also:
%   kpca, kernel_matrix, RBF_kernel


% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

if ~exist('fminunc'),
  error('This function needs the optimization function ''fminunc''.');
end    


if isstr(A1),
  kernel = A1;
  par = A2;
  eval('if isempty(Xt), Xt = A3; end','Xt = Xo;');
  eval('etype = A4;','etype = ''svd'';');
  eval('nb = A5;','nb = ''inf'';');
  [~,U] = kpca(Xo,kernel,par);
  [lam,U] = kpca(Xo,kernel,par,[], etype, nb);
else
  U = A1;
  lam = A2;
  kernel = A3;
  par = A4;
  eval('Xt = A5;','Xt = Xo;');
end

warning off
[nb,d] = size(Xt);
for n=1:nb,
  x = Xt(n,:);
  %dist_phi(x,x,kernel,par,U,lam,Xo)
  Xd(n,:) = fminunc(@dist_phi,x,[],x,kernel,par,U,lam,Xo);
end
warning on


function d = dist_phi(x,xor,kernel,par,U,lam,Xo)
% the distance in feature space between x and the subspace spanned
% by U,lam

k = kernel_matrix(Xo,kernel,par,[x;xor]);
betas = k(:,1)'*U;
alphas = k(:,2)'*U;
d = feval(kernel,x,x,par)-2*betas*alphas';