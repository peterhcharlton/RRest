function [cost,costs] = rcrossvalidate(model, L, wfun, estfct,combinefct)

% Estimate the model performance of a model with robust [$ l$] -fold crossvalidation

% CAUTION!! Use this function only to obtain the value of the rcrossvalidation score 
% function given the tuning parameters. Do not use this function together with 
% 'tunelssvm', but use 'rcrossvalidatelssvm' instead. The latter is a faster 
% implementation which uses previously computed results.
%
%
% >> cost = rcrossvalidate({Xtrain,Ytrain,type,gam,sig2})
% >> cost = rcrossvalidate( model)
%
% Robustness in the $l$-fold crossvalidation score function is obtained by 
% iteratively reweighting schemes.
%
% This routine is computational intensive.
%
%
% Some commonly used criteria are:
%
% >> cost = rcrossvalidate(model, 10, 'whuber', 'mae')
% >> cost = rcrossvalidate(model, 10, 'whampel', 'mae')
% >> cost = rcrossvalidate(model, 10, 'wlogistic', 'mae')
% >> cost = rcrossvalidate(model, 10, 'wmyriad', 'mae')
%
% Full syntax
%
%     1. Using LS-SVMlab with the functional interface:
%
% >> [cost, costs] = rcrossvalidate({X,Y,type,gam,sig2,kernel,preprocess}, L, wfun, estfct, combinefct)
%
%       Outputs
%         cost          : Cost estimation of the L-fold cross validation
%         costs(*)      : L x 1 vector with costs estimated on the L different folds
%       Inputs
%         X             : Training input data used for defining the LS-SVM and the preprocessing
%         Y             : Training output data used for defining the LS-SVM and the preprocessing
%         type          : 'function estimation' ('f') or 'classifier' ('c')
%         gam           : Regularization parameter
%         sig2          : Kernel parameter (bandwidth in the case of the 'RBF_kernel')
%         kernel(*)     : Kernel type (by default 'RBF_kernel')
%         preprocess(*) : 'preprocess'(*) or 'original'
%         L(*)          : Number of folds (by default 10)
%         wfun(*)       : weighting scheme (by default: whuber)
%         estfct(*)     : Function estimating the cost based on the residuals (by default mse)
%         combinefct(*) : Function combining the estimated costs on the different folds (by default mean)
%
%
%     2. Using the object oriented interface:
%
% >> [cost, costs] = crossvalidate(model, L, wfun, estfct, combinefct)
%
%       Outputs
%         cost          : Cost estimation of the L-fold cross validation
%         costs(*)      : L x 1 vector with costs estimated on the L different folds
%         ec(*)         : N x 1 vector with residuals of all data
%       Inputs
%         model         : Object oriented representation of the LS-SVM model
%         Xval          : Nt x d matrix with the inputs of the validation points used in the procedure
%         Yval          : Nt x m matrix with the outputs of the validation points used in the procedure
%         L(*)          : Number of folds (by default 10)
%         wfun(*)       : weighting scheme (by default: whuber)
%         estfct(*)     : Function estimating the cost based on the residuals (by default mse)
%         combinefct(*) : Function combining the estimated costs on the different folds (by default mean)
%
% See also:
% mae,whuber,wlogistic,whampel,wmyriad, crossvalidate, trainlssvm, robustlssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


% LS-SVMlab
eval('model = initlssvm(model{:});',' ');
eval('L;','L=min(ceil(sqrt(model.nb_data)),10);');
eval('estfct;','estfct=''mae'';');
eval('combinefct;','combinefct=''mean'';');
eval('wfun;','wfun=''whuber'';');
%
% initialisation and defaults
%
nb_data = size(model.ytrain,1);

if L==nb_data, p = 1:nb_data; else p = randperm(nb_data); end
px = model.xtrain(p,:);
py = model.ytrain(p,:);

[~,Y] = postlssvm(model,[],py);

%initialize: no incremental  memory allocation
costs = zeros(L,length(model.gam));
block_size = floor(nb_data/L);

S = ones(nb_data,1);
Atot = kernel_matrix(px,model.kernel_type,model.kernel_pars)+eye(nb_data)./model.gam;
%
%
% start loop over l validations
%
for l = 1:L,
    
    % divide in data and validation set, trainings data set is a copy
    % of permutated_data, validation set is just a logical index
    if l==L,
        train = 1:block_size*(l-1);
        validation = block_size*(l-1)+1:nb_data;
    else
        train = [1:block_size*(l-1) block_size*l+1:nb_data];
        validation = block_size*(l-1)+1:block_size*l;
    end
    
    A = [0 S(train)';S(train) Atot(train,train)];
    b = [0;py(train)];
    
    % Solve linear system
    sol = linsolve(A,b,struct('SYM',true));
    
    % Determine residuals ek
    ek = sol(2:end)./model.gam;
    g = model.gam;
    %for i=2:size(A,1), A(i,i) = A(i,i) - 1/g; end
    A = A-eye(size(train,2)+1)./g; A(1,1)=0;
    Ah = A;
    %
    % robust estimation of the variance
    %
    for k = 1:20
        vare = 1.483*median(abs((ek)-median(ek)));
        alphaold = sol(2:end);
        %
        % robust re-estimation of the alpha's and the b
        %
        cases = reshape((ek./vare),1,size(ek,1));
        W = g*weightingscheme(cases,wfun);
        
        for t=1:size(train,2), A(t+1,t+1) = A(t+1,t+1)+1./W(t); end
        
        sol = linsolve(A,b,struct('SYM',true));
        
        ek = sol(2:end)./W';
        A = Ah;
        if norm(abs(alphaold-sol(2:end)),'fro')<=1e-4,
            %fprintf('\n Converged after %.0f iteration(s)', k);
            k = inf;
        end
        model.status = 'changed';
    end
    
    % regression
    % Simulate system on validation data
    yh = Atot(train,validation)'*sol(2:end) + ones(numel(validation),1)*sol(1);
    [~,yh] = postlssvm(model,[],yh);
    costs(l,1) = feval(estfct,yh - Y(validation,:));
end
cost = feval(combinefct, costs);


