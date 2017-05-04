function [cost,costs] = crossvalidate(model, L, estfct,combinefct)

% Estimate the model performance of a model with [$ l$] -fold crossvalidation
%
% CAUTION!! Use this function only to obtain the value of the crossvalidation score
% function given the tuning parameters. Do not use this function together with
% 'tunelssvm', but use 'crossvalidatelssvm' instead. The latter is a faster
% implementation which uses previously computed results.
%
% >> cost = crossvalidate({Xtrain,Ytrain,type,gam,sig2})
% >> cost = crossvalidate( model)
%
% The data is once permutated randomly, then it is divided into L (by default 10)
% disjoint sets. In the i-th (i=1,...,l) iteration, the i-th set is used to estimate
% the performance ('validation set') of the model trained on the other l-1 sets ('training set').
% Finally, the l (denoted by L) different estimates of the performance are combined (by default by the 'mean').
% The assumption is made that the input data are distributed independent and identically over the
% input space. As additional output, the costs in the different folds ('costs') of the data are returned:
%
% >> [cost, costs] = crossvalidate(model)
%
% Some commonly used criteria are:
%
% >> cost = crossvalidate(model, 10, 'misclass', 'mean')
% >> cost = crossvalidate(model, 10, 'mse', 'mean')
% >> cost = crossvalidate(model, 10, 'mae', 'median')
%
% Full syntax
%
%     1. Using LS-SVMlab with the functional interface:
%
% >> [cost, costs] = crossvalidate({X,Y,type,gam,sig2,kernel,preprocess}, L, estfct, combinefct)
%
%       Outputs
%         cost          : Cost estimation of the L-fold cross validation
%         costs(*)      : L x 1 vector with costs estimated on the L different folds
%
%       Inputs
%         X             : Training input data used for defining the LS-SVM and the preprocessing
%         Y             : Training output data used for defining the LS-SVM and the preprocessing
%         type          : 'function estimation' ('f') or 'classifier' ('c')
%         gam           : Regularization parameter
%         sig2          : Kernel parameter (squared bandwidth in the case of the 'RBF_kernel')
%         kernel(*)     : Kernel type (by default 'RBF_kernel')
%         preprocess(*) : 'preprocess'(*) or 'original'
%         L(*)          : Number of folds (by default 10)
%         estfct(*)     : Function estimating the cost based on the residuals (by default mse)
%         combinefct(*) : Function combining the estimated costs on the different folds (by default mean)
%
%
%     2. Using the object oriented interface:
%
% >> [cost, costs] = crossvalidate(model, L, estfct, combinefct)
%
%       Outputs
%         cost          : Cost estimation of the L-fold cross validation
%         costs(*)      : L x 1 vector with costs estimated on the L different folds
%
%       Inputs
%         model         : Object oriented representation of the LS-SVM model
%         L(*)          : Number of folds (by default 10)
%         estfct(*)     : Function estimating the cost based on the residuals (by default mse)
%         combinefct(*) : Function combining the estimated costs on the different folds (by default mean)
%
%
%     3. Using other modeling techniques::
%
%
% See also:
% leaveoneout, gcrossvalidate, trainlssvm, simlssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

% LS-SVMlab
eval('model = initlssvm(model{:});',' ');
eval('L;','L=min(ceil(sqrt(model.nb_data)),10);');
eval('estfct;','estfct=''mse'';');
eval('combinefct;','combinefct=''mean'';');

%
% initialisation and defaults
%
nb_data = size(model.ytrain,1);
d = size(model.xtrain,2);

if L==nb_data, p = 1:nb_data; else p = randperm(nb_data); end
px = model.xtrain(p,:);
py = model.ytrain(p,:);

[~,Y] = postlssvm(model,[],py);

%initialize: no incremental  memory allocation
costs = zeros(L,length(model.gam));
block_size = floor(nb_data/L);

% calculate matrix for LS-SVM once for the entire data
S = ones(nb_data,1);
Inb = eye(nb_data);
K = kernel_matrix(px,model.kernel_type,model.kernel_pars);
Atot = K+Inb./model.gam;

% Cholesky factor
try R = chol(Atot);
    % Solve full system
    q = R\(R'\[py S]);
    p = q(:,2); q = q(:,1);
    s  = 1/sum(p);
    bias  = s*sum(q);
    alpha  = q - p*bias;
    
    % Two expensive steps yet more efficient that using LINSOLVE on each fold
    Ri = R\Inb;
    C = Ri*Ri' - s*(p)*p';
    
catch %R = cholinc(sparse(Atot),1e-5);
    A = [K+Inb./model.gam S; S' 0];
    C = pinv(A);
    alpha = C*[py;0];
    %bias = alpha(nb_data+1);
    alpha = alpha(1:nb_data);
end

% Solve full system
q = R\(R'\[py S]);
p = q(:,2); q = q(:,1);
s  = 1/sum(p);
bias  = s*sum(q);
alpha  = q - p*bias;

% Two expensive steps yet more efficient that using LINSOLVE on each fold
Ri = R\Inb;
C = Ri*Ri' - s*(p)*p';

% start loop over l validations
for l = 1:L,
    % divide data in validation set and trainings data set
    if l==L,
        validation = block_size*(l-1)+1:nb_data;
    else
        validation = block_size*(l-1)+1:block_size*l;
    end
    % Submatrix of C to compute residuals for the l-th fold left out
    Ckk = C(validation,validation);
    % Solution of small linear system (block_size x block_size)
    try % faster
        Rkk = chol(Ckk+eps);
        betak = Rkk\(Rkk'\alpha(validation));
    catch
        betak = Ckk\alpha(validation);
    end
    % latent outputs for validation
    yh = py(validation) - betak;
    [~,yh] = postlssvm(model,[],yh);
    if ~(model.type(1)=='c')
        costs(l,1) = feval(estfct,yh - Y(validation,:));
    else
        costs(l,1) = feval(estfct,Y(validation,:),sign(yh));
    end
end
cost = feval(combinefct,costs);