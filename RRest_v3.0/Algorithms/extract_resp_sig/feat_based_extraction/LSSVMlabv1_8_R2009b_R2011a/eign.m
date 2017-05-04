function [V,D,Ann] = eign(A1, A2,A3, args)
% Find the principal eigenvalues and eigenvectors of a matrix with Nystr�m's low rank approximation method
% 
% >>  D     = eign(A, nb)
% >> [V, D] = eign(A, nb)
% 
% In the case of using this method for low rank approximation and
% decomposing the kernel matrix, one can call the function without
% explicit construction of the matrix A.
% 
% >>  D     = eign(X, kernel, kernel_par, nb)
% >> [V, D] = eign(X, kernel, kernel_par, nb)
% 
%
% Full syntax
%  (We denote the size of positive definite matrix A with a*a.)
%
%     1. Given the full matrix:
% 
% >>  D    = eign(A,nb)
% >> [V,D] = eign(A,nb)
% 
%       
%       Outputs    
%         V(*)  : a x nb matrix with estimated principal eigenvectors of A
%         D     : nb x 1 vector with principal estimated eigenvalues of A
%       Inputs    
%         A     : a*a positive definite symmetric matrix
%         nb(*) : Number of approximated principal eigenvalues/eigenvectors
% 
%
%     2. Given the function to calculate the matrix elements:
% 
% >>  D = eign(X, kernel, kernel_par, nb)
% >> [V,D] = eign(X, kernel, kernel_par, nb)
% 
%       Outputs    
%         V(*)       : a x nb matrix with estimated principal eigenvectors of A
%         D          : nb x 1 vector with estimated principal eigenvalues of A
%       Inputs    
%         X          : N x d matrix with the training data
%         kernel     : Kernel type (e.g. 'RBF_kernel')
%         kernel_par : Kernel parameter (bandwidth in the case of the 'RBF_kernel')
%         nb(*)      : Number of eigenvalues/eigenvectors used in the eigenvalue decomposition approximation
% 
% See also:
%   eig, eigs, kpca, bay_lssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

% AFUN?
 if nargin~=2
  X = A1;
  kernel = A2;
  kernel_par = A3;
  N = size(X,1);
  
  %eval(['if args<1, error(''strict positive number of eigenvalues required;'');else n = args; end;'],...
  %'n=min(6,ceil(N*.75));');
    if args<1
      error('Strict positive number of eigenvalues required');
  else
      n = args;
  end;
  
  % random sampling
  s = randperm(N); sr=s(n+1:end); s=s(1:n); 
  %s = ceil(1:(N-1)/(n-1):N); s = s(1:n);
  
  ANn = zeros(n,N);
  ANn = kernel_matrix(X,kernel,kernel_par,X(s,:));
  Ann = ANn(s,:);

  % centering of matrix
  Zc = eye(n) - 1/n;
  %ZC = eye(N) - 1/N;
  Ann = Zc*Ann*Zc;
  %ANn = ZC*ANn*Zc;
  ANn = (Zc*(ANn*Zc)')';
else
  A = A1;
  N = size(A,1);
  
  %eval(['if args<1, error(''strict positive number of eigenvalues required;'');else n = args; end;'],...
  %'n=min(6,ceil(N*.75));');
  
  if A2<1
      error('Strict positive number of eigenvalues required');
  else
      n = A2;
  end;
  
  % random sampling
  s = randperm(N); sr=s(n+1:end); s=s(1:n); 
  %s = ceil(1:(N-1)/(n-1):N); s = s(1:n);

  
  ANn = A(:,s);
  Ann = A(s,s);
end




%
% compute eigenvalues en vectors of low rank approximation
%
[Vn,Dn] = eig(Ann);Dn = diag(Dn);


%
% select only relevant eigenvalues and sort
% (only largest eigenvectors are orthogonal)
%
[Dn,peff] = sort(Dn(find(Dn>1000*eps)));
Dn = Dn(end:-1:1); peff = peff(end:-1:1);
Vn = Vn(:,peff);

%
% Nystrom correction
%
D = (N/n).*Dn;


%
% eigenvectoren correctie
%
if nargout>1,
  %V = zeros(n,length(peff));
  for i=1:length(peff),
    V(:,i) = (sqrt(n/N)/Dn(peff(i)))*ANn*Vn(:,peff(i));    
  end
  
  %
  % correction of found eigenvectors:
  % orthogonal and unit length
  %

  % svd
  %[V,D2,ff] = svd(V*diag(D.^.5));   V = V(:,1:length(peff));
  
  % gram schmidt
  %[V,r] = gramschmidt2(V);
  %D = D.*r;

  %D = diag(D.^2);
else
  V=D;
end

