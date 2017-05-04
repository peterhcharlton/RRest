function cost = rcrossvalidatelssvm(model,Y, L, omega, estfct,combinefct)
%%%%%%%%%%%%%%%%%%%%%
% INTERNAL FUNCTION %
%%%%%%%%%%%%%%%%%%%%%
% Estimate the model performance of a model with l-fold robust crossvalidation

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @
% http://www.esat.kuleuven.be/sista/lssvmlab

%
% initialisation and defaults
%
%if size(X,1)~=size(Y,1), error('X and Y have different number of datapoints'); end
[nb_data,y_dim] = size(Y);
d = size(model.xtrain,2);
% LS-SVMlab
eval('model = initlssvm(model{:});',' ');
model.status = 'changed';

eval('L;','L=min(round(sqrt(size(model.xfull,1))),10);');
eval('estfct;','estfct=''mse'';');
eval('combinefct;','combinefct=''mean'';');

py = Y;
[~,Y] = postlssvm(model,[],Y);

gams = model.gamcsa; 
eval('sig2s = model.kernel_parscsa;','sig2s=[];')
eval('deltas = model.deltacsa;','deltas=[];')
%
%initialize: no incremental  memory allocation
%
costs = zeros(L,length(gams));
block_size = floor(nb_data/L);

% check whether there are more than one gamma or sigma
for j =1:numel(gams)
    if strcmp(model.kernel_type,'RBF_kernel') || strcmp(model.kernel_type,'RBF4_kernel')
        model = changelssvm(changelssvm(model,'gam',gams(j)),'kernel_pars',sig2s(j));
        eval('model.delta=deltas(j);','model.delta=[];')
    elseif strcmp(model.kernel_type,'lin_kernel')
        model = changelssvm(model,'gam',gams(j));
        eval('model.delta=deltas(j);','model.delta=[];')
    elseif strcmp(model.kernel_type,'poly_kernel')
        model = changelssvm(changelssvm(model,'gam',gams(j)),'kernel_pars',[sig2s(1,j);sig2s(2,j)]);
        eval('model.delta=deltas(j);','model.delta=[];')
    else
        model = changelssvm(changelssvm(model,'gam',gams(j)),'kernel_pars',[sig2s(1,j);sig2s(2,j);sig2s(3,j)]);
        eval('model.delta=deltas(j);','model.delta=[];')
    end
    
    S = ones(nb_data,1); 
    Atot = kernel_matrix2(omega,model.kernel_type,model.kernel_pars,d)+eye(nb_data)./model.gam;
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
        for k = 1:30
            vare = 1.483*median(abs((ek)-median(ek)));
            alphaold = sol(2:end);
            %
            % robust re-estimation of the alpha's and the b
            %
            cases = reshape((ek./vare),1,size(ek,1));
            W = g*weightingscheme(cases,model.weights,model.delta);
            
            for t=1:size(train,2), A(t+1,t+1) = A(t+1,t+1)+1./W(t); end   
                        
            sol = linsolve(A,b,struct('SYM',true));
                        
            ek = sol(2:end)./W';
            A = Ah;
            if norm(abs(alphaold-sol(2:end)),'fro')<=1e-3,
                %fprintf('\n Converged after %.0f iteration(s)', k);
                k = inf;                
            end
            model.status = 'changed';
        end
        
        % regression
        % Simulate system on validation data
        yh = Atot(train,validation)'*sol(2:end) + ones(numel(validation),1)*sol(1);
        [~,yh] = postlssvm(model,[],yh);
        z = yh - Y(validation,:);
        eval('costs(l,j) = feval(estfct,z);')
    end
end

cost = feval(combinefct, costs);
        
        
    