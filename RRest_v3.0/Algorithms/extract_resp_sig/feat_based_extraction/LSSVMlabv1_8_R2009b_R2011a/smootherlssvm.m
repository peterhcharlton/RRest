function S = smootherlssvm(model,Xt)

% Calculates the smoother matrix for LS-SVM.
% Inputs:
%        - model : object oriented interface of the FS-LSSVM model
%        - Xt (*): smoother matrix evaluated in test point(s) Xt
%
% Outputs:
%        - S: smoother matrix

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


if isempty(model.gam) && isempty(model.kernel_pars)
    error('Please supply one or more learning parameters');
end

K = kernel_matrix(model.xtrain,model.kernel_type,model.kernel_pars);
if nargin < 2
    S = smoother(model,K);
else
    if model.preprocess(1)=='p'
        % Preprocess the test data
        Xt = prelssvm(model,Xt);
    end
    Kt = kernel_matrix(model.xtrain,model.kernel_type,model.kernel_pars,Xt);
    S = smoother(model,K,Kt);
end


function S = smoother(model,K,varargin)
Z = pinv(K+eye(model.nb_data)./model.gam);
c = sum(sum(Z));
J = ones(model.nb_data)./c;
if isempty(varargin)
    S = K*(Z-Z*J*Z) + J*Z;
else
    Kt = varargin{1}';
    J1 = ones(size(Kt,1),size(Z,1))./c;
    S = Kt*(Z-Z*J*Z) + J1*Z;
end
