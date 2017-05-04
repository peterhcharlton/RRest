function [eigval, eigvec, scores, omega,recErrors,optOut] = kpca(Xtrain, kernel_type, kernel_pars ,Xt,etype,nb,rescaling)
% Kernel Principal Component Analysis (KPCA)
%
% >> [eigval, eigvec] = kpca(X, kernel_fct, sig2)
% >> [eigval, eigvec, scores] = kpca(X, kernel_fct, sig2, Xt)
%
% Compute the nb largest eigenvalues and the corresponding rescaled
% eigenvectors corresponding with the principal components in the
% feature space of the centered kernel matrix. To calculate the
% eigenvalue decomposition of this N x N matrix, Matlab's
% eig is called by default. The decomposition can also be
% approximated by Matlab ('eigs') or by Nystrom method ('eign')
% using nb components. In some cases one wants to disable
% ('original') the rescaling of the principal components in feature
% space to unit length.
%
% The scores of a test set Xt on the principal components is computed by the call:
%
% >> [eigval, eigvec, scores] = kpca(X, kernel_fct, sig2, Xt)
%
% Full syntax
%
% >> [eigval, eigvec, empty, omega] = kpca(X, kernel_fct, sig2)
% >> [eigval, eigvec, empty, omega] = kpca(X, kernel_fct, sig2, [],etype, nb)
% >> [eigval, eigvec, empty, omega] = kpca(X, kernel_fct, sig2, [],etype, nb, rescaling)
% >> [eigval, eigvec, scores, omega] = kpca(X, kernel_fct, sig2, Xt)
% >> [eigval, eigvec, scores, omega] = kpca(X, kernel_fct, sig2, Xt,etype, nb)
% >> [eigval, eigvec, scores, omega] = kpca(X, kernel_fct, sig2, Xt,etype, nb, rescaling)
% >> [eigval, eigvec, scores, omega, recErrors] = kpca(X, kernel_fct, sig2, Xt)
% >> [eigval, eigvec, scores, omega, recErrors] = kpca(X, kernel_fct, sig2, Xt,etype, nb)
% >> [eigval, eigvec, scores, omega, recErrors] = kpca(X, kernel_fct, sig2, Xt,etype, nb, rescaling)
% >> [eigval, eigvec, scores, omega, recErrors, optOut] = kpca(X, kernel_fct, sig2, Xt)
% >> [eigval, eigvec, scores, omega, recErrors, optOut] = kpca(X, kernel_fct, sig2, Xt,etype, nb)
% >> [eigval, eigvec, scores, omega, recErrors, optOut] = kpca(X, kernel_fct, sig2, Xt,etype, nb, rescaling)
%
% Outputs
%   eigval       : N (nb) x 1 vector with eigenvalues values
%   eigvec       : N x N (N x nb) matrix with the principal directions
%   scores(*)    : Nt x nb matrix with the scores of the test data (or [])
%   omega(*)     : N x N centered kernel matrix
%   recErrors(*) : Nt x 1 vector with the reconstruction error of the test data
%   optOut(*)    : Optional cell array containing the centered test kernel matrix in optOut{1}
%                  and the squared 2-norm of the test points in the feature space in optOut{2}
% Inputs
%   X            : N x d matrix with the inputs of the training data
%   kernel       : Kernel type (e.g. 'RBF_kernel')
%   sig2         : Kernel parameter(s) (for linear kernel, use [])
%   Xt(*)        : Nt x d matrix with the inputs of the test data (or [])
%   etype(*)     : 'svd', 'eig'(*),'eigs','eign'
%   nb(*)        : Number of eigenvalues/eigenvectors used in the eigenvalue decomposition approximation
%   rescaling(*) : 'original size' ('o') or 'rescaled'(*) ('r')
%
% See also:
%   bay_lssvm, bay_optimize, eign

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


%
% defaults
%
nb_data = size(Xtrain,1);


if ~exist('nb','var')
    nb=10;
end;

if ~exist('etype','var')
    etype='eig';
end;

if ~exist('rescaling','var')
    rescaling='r';
end;




%eval('n=min(n,nb_data);','n=min(10,nb_data);')
% eval('centr;','centr=''rescaled'';');
% eval('etype;','etype=''eig'';');
% eval('Xt;','Xt=[];');


%
% tests
%
if exist('Xt','var') && ~isempty(Xt) && size(Xt,2)~=size(Xtrain,2),
    error('Training points and test points need to have the same dimension');
end

if ~(strcmpi(etype,'svd') || strcmpi(etype,'eig') || strcmpi(etype,'eigs') || strcmpi(etype,'eign')),
    error('Eigenvalue decomposition via ''svd'', ''eig'', ''eigs'' or ''eign''...');
end


if (strcmpi(etype,'svd') || strcmpi(etype,'eig') || strcmpi(etype,'eigs')),
    
    omega = kernel_matrix(Xtrain,kernel_type, kernel_pars);
    
    % Centering
    Meanvec = mean(omega,2);
    MM = mean(Meanvec);
    omega=omega-Meanvec*ones(1,nb_data)-ones(nb_data,1)*Meanvec'+MM;
    
    % numerical stability issues
    %omega = (omega+omega')./2;
    
    if strcmpi(etype,'svd'),
        [eigvec, eigval] = svd(omega);
    elseif strcmpi(etype,'eig'),
        [eigvec, eigval] = eig(omega);
    elseif (strcmpi(etype,'eigs')),
        [eigvec, eigval] = eigs(omega,nb);
    end
    %eigval = diag(eigval)./(nb_data-1);
    eigval=diag(eigval);
    
    
    
elseif strcmpi(etype,'eign'),
    if nargout>1,
        [eigvec,eigval] = eign(Xtrain,kernel_type,kernel_pars, nb);
    else
        eigval = eign(Xtrain,kernel_type,kernel_pars, nb);
    end
    omega = [];
    %eigval = (eigval)./(nb_data-1);
    Meanvec = [];
    MM = [];
    
else
    error('Unknown type for eigenvalue approximation');
end


%% Eigenvalue/vector sorting in descending order

[eigval,evidx]=sort(eigval,'descend');
eigvec=eigvec(:,evidx);



%%


%
% only keep relevant eigvals & eigvec
%
peff = find(eigval>1000*eps);
%eigval = eigval(peff);
neff = length(peff);
%if nargout>1, eigvec = eigvec(:,peff); end

% rescaling the eigenvectors
if (rescaling(1) =='r' && nargout>1),
    %disp('rescaling the eigvec');
    for i=1:neff,
        eigvec(:,i) = eigvec(:,i)./sqrt(eigval(i));
    end
end


%
% compute scores
%
if exist('Xt','var') && ~isempty(Xt),
    nt=size(Xt,1);
    omega_t = kernel_matrix(Xtrain,kernel_type, kernel_pars,Xt);
    MeanvecT=mean(omega_t,1);
    omega_t=omega_t-Meanvec*ones(1,nt) - ones(nb_data,1)*MeanvecT+MM;
    scores = omega_t'*eigvec;
    
    normProjXt=diag(omega_t'*(eigvec*eigvec')*omega_t);
    
    
    if strcmp(kernel_type,'RBF_kernel')
        ks = ones(1,nt);
    else
        
        for i = 1:nt
            ks(i) = feval(kernel_type,Xt(i,:),Xt(i,:),kernel_pars);
        end;
        
    end;
    
    normPhiXt=ks'-MeanvecT'*2+MM;
    
    recErrors= normPhiXt-normProjXt;
    
    optOut={omega_t,normPhiXt};
else
    scores = [];
end




