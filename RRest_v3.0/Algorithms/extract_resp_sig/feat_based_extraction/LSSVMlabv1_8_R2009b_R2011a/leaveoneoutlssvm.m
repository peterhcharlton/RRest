function cost = leaveoneoutlssvm(model,Y,omega, estfct)
% Fast leave-one-out cross-validation for the LS-SVM based on one full matrix inversion
%
%%%%%%%%%%%%%%%%%%%%%
% INTERNAL FUNCTION %
%%%%%%%%%%%%%%%%%%%%%
% Estimate the model performance of a model with fast LOO crossvalidation.
% This implementation is based on one full matrix inverse. Implementation
% based on "Z. Ying and K.C. Keong: Fast Leave-One-Out Evaluation and
% Improvement on Inference for LS-SVM's, Proc. ICPR, 2004"

% Copyright (c) 2010,  KULeuven-ESAT-SCD, License & help @% http://www.esat.kuleuven.ac.be/sista/lssvmlab

%
% See also:
%   leaveoneout, crossvalidate, trainlssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

% LS-SVMlab
eval('model = initlssvm(model{:});',' ');
model.status = 'changed';
d = size(model.xtrain,2);

eval('estfct;','estfct=''mse'';');
eval('combinefct;','combinefct=''mean'';');

gams = model.gamcsa; try sig2s = model.kernel_parscsa; catch, sig2s = [];end

%
%initialize: no incremental  memory allocation
%
cost = zeros(1,length(gams));
py = Y;
[~,Y] = postlssvm(model,[],Y); % Y is raw data, non preprocessed

% check whether there are more than one gamma or sigma
for g =1:numel(gams)
    if strcmp(model.kernel_type,'RBF_kernel') || strcmp(model.kernel_type,'RBF4_kernel')
        model = changelssvm(changelssvm(model,'gam',gams(g)),'kernel_pars',sig2s(g));
    elseif strcmp(model.kernel_type,'lin_kernel')
        model = changelssvm(model,'gam',gams(g));
    elseif strcmp(model.kernel_type,'poly_kernel')
        model = changelssvm(changelssvm(model,'gam',gams(g)),'kernel_pars',[sig2s(1,g);sig2s(2,g)]);
    else
        model = changelssvm(changelssvm(model,'gam',gams(g)),'kernel_pars',[sig2s(1,g);sig2s(2,g);sig2s(3,g)]);
    end
    
    % kernel matrix computation
    K = kernel_matrix2(omega,model.kernel_type,model.kernel_pars,d);
    
    Ka = pinv([K+eye(model.nb_data)./model.gam ones(model.nb_data,1);ones(1,model.nb_data) 0]);
    sol = Ka*[py;0]; model.alpha = sol(1:end-1); model.b = sol(end);
    yh = py - model.alpha./diag(Ka(1:model.nb_data,1:model.nb_data));
    
    [~,yh] = postlssvm(model,[],yh);
    if ~(model.type(1)=='c')
        cost(g) = feval(estfct,yh-Y);
    else
        cost(g) = feval(estfct,Y,sign(yh));
    end
    
end
