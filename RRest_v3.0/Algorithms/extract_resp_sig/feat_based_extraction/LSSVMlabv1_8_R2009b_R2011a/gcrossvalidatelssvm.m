function cost = gcrossvalidatelssvm(model,Y,omega,estfct)

%%%%%%%%%%%%%%%%%%%%%
% INTERNAL FUNCTION %
%%%%%%%%%%%%%%%%%%%%%
% Estimate the model performance of a model with genralized crossvalidation
% for regression with the LS-SVM

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @
% http://www.esat.kuleuven.be/sista/lssvmlab

% LS-SVMlab
eval('model = initlssvm(model{:});',' ');
model.status = 'changed';
d = size(model.xtrain,2);

eval('estfct;','estfct=''mse'';');

gams = model.gamcsa; try sig2s = model.kernel_parscsa; catch, sig2s = [];end

py = Y;
[~,Y] = postlssvm(model,[],Y); % Y is raw data, non preprocessed

cost = zeros(length(gams),1);
% check whether there are more than one gamma or sigma
for j =1:numel(gams)
    if strcmp(model.kernel_type,'RBF_kernel') || strcmp(model.kernel_type,'RBF4_kernel')
        model = changelssvm(changelssvm(model,'gam',gams(j)),'kernel_pars',sig2s(j));
    elseif strcmp(model.kernel_type,'lin_kernel')
        model = changelssvm(model,'gam',gams(j));
    else
        model = changelssvm(changelssvm(model,'gam',gams),'kernel_pars',[sig2s(1,j);sig2s(2,j)]);
    end
    
    % Calculate Kernel matrix and trace smoother
    K = kernel_matrix2(omega,model.kernel_type,model.kernel_pars,d);
    tr = trace(smoother(model,K));
    
    % Solve the linear system
    S = ones(model.nb_data,1);
    sol = linsolve([0 S';S K+eye(model.nb_data)./model.gam],[0;py],struct('SYM',true));
    
    % Simulate
    %ek = sol(2:end)./model.gam;
    yh = K*sol(2:end) + ones(model.nb_data,1)*sol(1);
    [~,yh] = postlssvm(model,[],yh);
    % Generalized cross-validation
    if ~(model.type(1)=='c')
        cost(j,1) = feval(estfct,yh-Y);
    else
        cost(j,1) = feval(estfct,Y,sign(yh));
    end
    cost(j,1) = cost(j,1)/((1-tr/size(model.ytrain,1))^2);
end

function S = smoother(model,K)
% Smoother Matrix of the LS-SVM
%
% f = K*alpha+1*b
% f = S*y
Z = pinv(K+eye(model.nb_data)./model.gam);
c = sum(sum(Z));
J = ones(model.nb_data)./c;
S = K*(Z-Z*J*Z) + J*Z;



