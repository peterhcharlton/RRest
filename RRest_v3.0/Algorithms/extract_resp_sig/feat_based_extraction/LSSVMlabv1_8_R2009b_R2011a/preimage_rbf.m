function Ximg = preimage_rbf(Xtr,sig2,U,B,type,npcs,maxIts)
%
% function Ximg = preimage_rbf(Xtr,sig2,U,B,type,npcs,maxIts)
%   Reconstruction or denoising after kernel PCA with RBF kernels, i.e. to find the
%   approximate preimage (in the input space) of the corresponding feature space expansions
% 
% Inputs 
%   Xtr     : N by d matrix of the training data used to find the prinicipal components.
%   sig2    : parameter for the RBF kernel used, k(x,z)=exp(-norm(x-z)^2/sig2).
%   U       : the eigenvectors computed from the kernel PCA using RBF kernel with parameter sig2.
%   B       : for reconstruction, B is the compressed data,
%               i.e. the projections of the data on to the first n PCs;
%             for denoising, B is the Nt by d matrix of original noisy data;
%             if not specified, Xtr is denoised instead.
%   type    : 'reconstruct' or 'denoise'
%   npcs    : number of PCs used for approximation
%   maxIts  : maximum iterations allowed to update the preimage, 1000 by default.
%
% Outputs
%   Ximg    : the reconstructed or denoised data in the input space
% 
% Usage e.g. 
%   >> [lam,U] = kpca(Xtr,'RBF_kernel',sig2);
%   >> [lam, perm] = sort(-lam); lam = -lam; U = U(:,perm); 
%   >> projections = kernel_matrix(Xtr,'RBF_kernel',sig2,Xtest)'*U;
%   >> Xr = preimage_rbf(Xtr,sig2,U,projections(:,1:npcs),'r'); % Reconstruction
%   >> Xd = preimage_rbf(Xtr,sig2,U(:,1:npcs),Xnoisy,'d');      % Denoising
%   >> Xdtr = preimage_rbf(Xtr,sig2,U(:,1:npcs));  % Denoising on the training data
%
% see also:
%    kpca, denoise_kpca, RBF_kernel
%
% Reference
%   Mika S., Schoelkopf B., Smola A., Muller K.-R., Scholz M., Ratsch G. (1999), ``Kernel
%   PCA and de-noising in feature spaces'', Advances in Neural Information Processing
%   Systems 11, 536-542, MIT Press. 
%

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

MAXDX=1e-6; % Convergence criterion
[~, dim]=size(Xtr);
mX=mean(Xtr); sX=std(Xtr);

if nargin<4, B = Xtr; end
[Nt, dimB] = size(B);

if nargin<5
    if dimB==dim, type='denoise'; else type='reconstruct'; end
else
    if type(1)~='d'&type(1)~='r'| (type(1)=='d'&dimB~=dim),
        warning('Invalid type specified, default value is used!')
        if dimB==dim, type='denoise'; else type='reconstruct'; end 
    end
end

if nargin<6
    if type(1)=='r'; npcs=dimB; else npcs=size(U,2); end
else
    if npcs>size(U,2) | (type(1)=='r' & npcs~=dimB),
        warning('Invalid number of PCs, default value is used!'),
        if type(1)=='r'; npcs=dimB; else npcs=size(U,2); end
    end
end

if nargin<7, maxIts=1000; end

U=U(:,1:npcs); 

for n=1:Nt
    cont=1; t=0; ts=0; 
    if type(1)=='r'
        % reconstuction
        rs = U*B(n,:)'; 
        x = zeros(1, dim); % set the initial value of approximate preimage for reconstruction
        k = RBF_kernel(x,Xtr,sig2);         
    else
        % denoise
        x = B(n,:); % initial value of the approximate preimage for denoising is the noisy data
        k = RBF_kernel(x,Xtr,sig2); 
        rs = U*(k'*U)';
    end                             

    % iteratively update the approximate preimage x
    %
    while cont, 

        d=rs'*k; % the reconstruction error is (const-2d)

        if d==0; 
            % choose a different starting value  
%            [k,id]=min(k); x_new = Xtr(id,:);
            randn('state',cputime+ts); x_new=(randn(1,dim)+mX).*sX;  
            fprintf('%5d> Starting value changed!(d=0) \n', ts);  
        else 
            % update approximate preimage with a linear combination of kpca training data
            x_new = sum(rs.*k*ones(1,dim).*Xtr)/d;
        end

        dx = norm(x_new - x); 
        x = x_new;

        t=t+1; ts=ts+1;
    
        if dx<MAXDX; 
            cont = 0; 
         % fprintf('%5d> Converged! \n', ts); 
        else
            if ts>=maxIts;
                cont=0;
                fprintf('%5d> Maximum iteration reached!\n', ts); 
            elseif t>=500;
                % choose a different starting value                 
%                [k,id]=min(k); x = Xtr(id,:); 
                randn('state',cputime+ts); x=(randn(1,dim)+mX).*sX;  
                t=0;
            end
        end
        k = RBF_kernel(x,Xtr,sig2); 
        
    end % while
    Ximg(n,:)=x_new;
end % for
