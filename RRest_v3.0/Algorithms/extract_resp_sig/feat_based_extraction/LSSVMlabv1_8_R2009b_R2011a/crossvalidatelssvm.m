function cost = crossvalidatelssvm(model,Y, L, omega, estfct,combinefct)
% Estimate the model performance of a model with l-fold crossvalidation
%
%%%%%%%%%%%%%%%%%%%%%
% INTERNAL FUNCTION %
%%%%%%%%%%%%%%%%%%%%%
% Estimate the model performance of a model with fast l-fold crossvalidation.
% Implementation based on "De Brabanter et al., Computationsl Statistics & Data Analysis, 2010"

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @% http://www.esat.kuleuven.be/sista/lssvmlab

%
% See also:
%   leaveoneoutlssvm, crossvalidatelssvm, trainlssvm
% Copyright (c) 2002,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.ac.be/sista/lssvmlab
% initialisation and defaults
%
%if size(X,1)~=size(Y,1), error('X and Y have different number of datapoints'); end
nb_data = size(Y,1);
d = size(model.xtrain,2);
% LS-SVMlab
eval('model = initlssvm(model{:});',' ');
model.status = 'changed';

eval('L;','L=min(round(sqrt(size(model.xfull,1))),10);');
eval('estfct;','estfct=''mse'';');
eval('combinefct;','combinefct=''mean'';');

% Y is raw data, non preprocessed
py = Y;
[~,Y] = postlssvm(model,[],Y);

gams = model.gamcsa; try sig2s = model.kernel_parscsa; catch, sig2s = [];end

%initialize: no incremental  memory allocation
costs = zeros(L,length(gams));
block_size = floor(nb_data/L);

% check whether there are more than one gamma or sigma
for j =1:numel(gams)
    if strcmp(model.kernel_type,'RBF_kernel') || strcmp(model.kernel_type,'RBF4_kernel')
        model = changelssvm(changelssvm(model,'gam',gams(j)),'kernel_pars',sig2s(j));
    elseif strcmp(model.kernel_type,'lin_kernel')
        model = changelssvm(model,'gam',gams(j));
    elseif strcmp(model.kernel_type,'poly_kernel')
        model = changelssvm(changelssvm(model,'gam',gams(j)),'kernel_pars',[sig2s(1,j);sig2s(2,j)]);
    else
        model = changelssvm(changelssvm(model,'gam',gams(j)),'kernel_pars',[sig2s(1,j);sig2s(2,j);sig2s(3,j)]);
    end
    
    
    % calculate matrix for LS-SVM once for the entire data
    S = ones(nb_data,1);
    Inb = eye(nb_data);
    K = kernel_matrix2(omega,model.kernel_type,model.kernel_pars,d);
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
    
    % start loop over l validations
    for l = 1:L,
        % divide data in validation set and trainings data set
        if l==L,
            %%train = 1:block_size*(l-1); % not used
            validation = block_size*(l-1)+1:nb_data;
        else
            %%train = [1:block_size*(l-1) block_size*l+1:nb_data]; % not used
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
            costs(l,j) = feval(estfct,yh - Y(validation,:));
        else
            costs(l,j) = feval(estfct,Y(validation,:),sign(yh));
        end
    end
end
cost = feval(combinefct, costs);