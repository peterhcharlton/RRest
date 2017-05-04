function [Pplus, Pmin, bay,model] = bay_modoutClass(model,X,priorpos,type,nb,bay)
% Estimate the posterior class probabilities of a binary classifier using Bayesian inference
%
% >> [Ppos, Pneg] = bay_modoutClass({X,Y,'classifier',gam,sig2}, Xt)
% >> [Ppos, Pneg] = bay_modoutClass(model, Xt)
% 
% Calculate the probability that a point will belong to the
% positive or negative classes taking into account the uncertainty
% of the parameters. Optionally, one can express prior knowledge as
% a probability between 0 and 1, where prior equal to 2/3 means
% that the  prior positive class probability is 2/3 (more likely to
% occur than the negative class).
% For binary classification tasks with a 2 dimensional input space,
% one can make a surface plot by replacing Xt by the string 'figure'.
% 
% Full syntax
% 
%     1. Using the functional interface:
% 
% >> [Ppos, Pneg] = bay_modoutClass({X,Y,'classifier',gam,sig2,kernel, preprocess}, Xt)
% >> [Ppos, Pneg] = bay_modoutClass({X,Y,'classifier',gam,sig2,kernel, preprocess}, Xt, prior)
% >> [Ppos, Pneg] = bay_modoutClass({X,Y,'classifier',gam,sig2,kernel, preprocess}, Xt, prior, type)
% >> [Ppos, Pneg] = bay_modoutClass({X,Y,'classifier',gam,sig2,kernel, preprocess}, Xt, prior, type, nb)
% >> bay_modoutClass({X,Y,'classifier',gam,sig2, kernel, preprocess}, 'figure')
% >> bay_modoutClass({X,Y,'classifier',gam,sig2, kernel, preprocess}, 'figure', prior)
% >> bay_modoutClass({X,Y,'classifier',gam,sig2, kernel, preprocess}, 'figure', prior, type)
% >> bay_modoutClass({X,Y,'classifier',gam,sig2, kernel, preprocess}, 'figure', prior, type, nb)
% 
%       Outputs    
%         Ppos    : Nt x 1 vector with probabilities that testdata Xt belong to the positive class
%         Pneg    : Nt x 1 vector with probabilities that testdata Xt belong to the negative(zero) class
%       Inputs    
%         X        : N x d matrix with the inputs of the training data
%         Y        : N x 1 vector with the outputs of the training data
%         type     : 'function estimation' ('f') or 'classifier' ('c')
%         gam      : Regularization parameter
%         sig2     : Kernel parameter (bandwidth in the case of the 'RBF_kernel')
%         kernel(*) : Kernel type (by default 'RBF_kernel')
%         preprocess(*) : 'preprocess'(*) or 'original'
%         Xt(*)    : Nt x d matrix with the inputs of the test data
%         prior(*) : Prior knowledge of the balancing of the training data (or [])
%         type(*)  : 'svd'(*), 'eig', 'eigs' or 'eign'
%         nb(*)    : Number of eigenvalues/eigenvectors used in the eigenvalue decomposition approximation
%
%     2. Using the object oriented interface:
% 
% >> [Ppos, Pneg, bay, model] = bay_modoutClass(model, Xt)
% >> [Ppos, Pneg, bay, model] = bay_modoutClass(model, Xt, prior)
% >> [Ppos, Pneg, bay, model] = bay_modoutClass(model, Xt, prior, type)
% >> [Ppos, Pneg, bay, model] = bay_modoutClass(model, Xt, prior, type, nb)
% >> bay_modoutClass(model, 'figure')
% >> bay_modoutClass(model, 'figure', prior)
% >> bay_modoutClass(model, 'figure', prior, type)
% >> bay_modoutClass(model, 'figure', prior, type, nb)
% 
%       Outputs    
%         Ppos     : Nt x 1 vector with probabilities that testdata Xt belong to the positive class
%         Pneg     : Nt x 1 vector with probabilities that testdata Xt belong to the negative(zero) class
%         bay(*)   : Object oriented representation of the results of the Bayesian inference
%         model(*) : Object oriented representation of the LS-SVM model
%       Inputs    
%         model    : Object oriented representation of the LS-SVM model
%         Xt(*)    : Nt x d matrix with the inputs of the test data
%         prior(*) :Prior knowledge of the balancing of the training data (or [])
%         type(*)  : 'svd'(*), 'eig', 'eigs' or 'eign'
%         nb(*)    : Number of eigenvalues/eigenvectors used in the eigenvalue decomposition approximation
% 
% See also:
%   bay_lssvm, bay_optimize, bay_errorbar, ROC


% Copyright (c) 2002,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.ac.be/sista/lssvmlab


% default handling
if iscell(model),
  model = trainlssvm(model);
end

if (model.type(1)~='c'), 
  error('this moderated output only possible for classification...'); 
end
eval('type;','type=''svd'';');
eval('nb;','nb=model.nb_data;');
if ~(strcmpi(type,'svd') | strcmpi(type,'eig') | strcmpi(type,'eigs') | strcmpi(type,'eign')),
  error('Eigenvalue decomposition via ''svd'', ''eig'', ''eigs'' or ''eign''...');
end
if strcmpi(type,'eign')
  warning('The resulting errorbars are most probably not very usefull...');  
end
eval('priorpos;','priorpos = .5*ones(model.y_dim,1);');
if isempty(priorpos), priorpos = .5*ones(model.y_dim,1); end
if ~isstr(X) & size(X,2)~=model.x_dim, 
  error('dimension datapoints is not equal to dimension of trainingspoints...');
end



if ~isstr(X),
    
  eval('[Pplus, Pmin, bay] = bay_modoutClassIn(model,X,priorpos,type,nb,bay);',...
       '[Pplus, Pmin, bay] = bay_modoutClassIn(model,X,priorpos,type,nb);');
  


% plot the curve including error bars
else
  if (model.x_dim==2 & model.y_dim==1),

    grain = 25;

    Xr = postlssvm(model,model.xtrain);
    disp(' COMPUTING PLOT OF MODERATED OUTPUT');

    % make grid
    Xmin = min(Xr,[],1);
    Xmax = max(Xr,[],1);
    Xs1 = (Xmin(1)):((Xmax(1)-Xmin(1))/grain):(Xmax(1));
    Xs2 = (Xmin(2)):((Xmax(2)-Xmin(2))/grain):(Xmax(2));
    grain = length(Xs1);
    
    [XX,YY] = meshgrid(Xs1,Xs2);
    l = size(XX,1)*size(XX,2);
    X = [reshape(XX,l,1) reshape(YY,l,1)];

    % compute moderated output 
    eval('[Pplus, Pmin, bay] = bay_modoutClassIn(model,X,priorpos,type,nb,bay);',...
	 '[Pplus, Pmin, bay] = bay_modoutClassIn(model,X,priorpos,type,nb);');

    
    figure;
    hold on;
    if isempty(model.kernel_pars),      
      title(['LS-SVM_{\gamma=' num2str(model.gam(1)) ...
             '}^{' model.kernel_type(1:3) '}, with moderated output' ...
             ' P_{pos} indicated by surface plot']);
    else
      title(['LS-SVM_{\gamma=' num2str(model.gam(1)) ', \sigma^2=' num2str(model.kernel_pars(1)) ...
             '}^{' model.kernel_type(1:3) '}, with moderated output' ...
             ' P_{pos} indicated by surface plot']);
    end
    xlabel('X_1');
    ylabel('X_2');
    zlabel('Y');
    surf(Xs1,Xs2,reshape(Pplus,grain,grain));
    
    
    % plot datapoints
    s = find(model.ytrain(:,1)>0);
    pp = plot3(Xr(s,1),Xr(s,2),ones(length(s),1) ,'*k');
    s = find(model.ytrain(:,1)<=0);
    pn = plot3(Xr(s,1),Xr(s,2),ones(length(s),1) ,'sk');
    legend([pp pn],'positive class','negative class');
    shading interp;
    colormap cool;
    axis([Xmin(1) Xmax(1) Xmin(2) Xmax(2)]);
    %colorbar
  else
    error(['cannot make a plot, give points to estimate confidence bounds instead...']);
  end
end






function [Pplus, Pmin, bay] = bay_modoutClassIn(model,X,priorpos,type, nb, bay)

% multiclass moderated output: recursive calls
if (model.y_dim>1), 
  %error('moderated output only possible for single class...'); 
  for i=1:model.y_dim,
    mff = model;
    mff.y_dim=1; 
    mff.ytrain=model.ytrain(:,i);
    mff.alpha = model.alpha(:,i);
    mff.b = model.b(i);
    mff.code='original';
    mff.preprocess = 'original';
    [Py(:,i), Pplus(:,i), Pmin(:,i), bay{i}] = bay_modoutClass(mff,X,priorpos(i),type,nb);
  end
  return
end


%
% evaluate LS-SVM in trainpoints, latent variables
%
Psv = latentlssvm(model,postlssvm(model,model.xtrain));
eval('Pymp = mean(Psv(find(Psv>0))));','Pymp=1;');
eval('Pymn = mean(Psv(find(Psv<=0)));','Pymp=-1;');

%model.latent  = 'no';
Py = latentlssvm(model,X);
nD = size(X,1);

% previous inference
eval('[FF1, FF2, FF3, bay] = bay_lssvm(model,1,type,nb);');

% kernel matrices
omega = kernel_matrix(model.xtrain,model.kernel_type, model.kernel_pars);
theta = kernel_matrix(model.xtrain,model.kernel_type, model.kernel_pars,X);
oo = ones(1,model.nb_data)*omega;
Zc = eye(model.nb_data) - ones(model.nb_data,1)*ones(1,model.nb_data)./model.nb_data;
Diagmatrix = (1/bay.mu - 1./(bay.zeta*bay.eigvals+bay.mu));

for i=1:nD,
  kxx(i,1) = feval(model.kernel_type, X(i,:),X(i,:), model.kernel_pars);
end

% positive class
  Mplusindex = (model.ytrain(:,1)>0);
  Nplus = sum(Mplusindex);
  Oplus = omega(:,Mplusindex);
  Oplusplus = omega(Mplusindex, Mplusindex);
  thetaplus = theta(Mplusindex,:);
  
  for i =1:nD,
    thetapluse(i,:) = (theta(:,i) - (1/Nplus)*sum(Oplus,2))'*Zc*bay.Rscores;
  end
  
  term1 = kxx - 2/(Nplus)*sum(thetaplus,1)';
  term2 = Nplus^-2 *sum(sum(Oplusplus));
  term3 = thetapluse.^2 * Diagmatrix;
  var_plus = (term1 + term2)./bay.mu - term3;
  
  
% negative class
  Mminindex = model.ytrain(:,1)<=0;
  Nmin = sum(Mminindex);
  Omin = omega(:,Mminindex);
  Ominmin = omega(Mminindex, Mminindex);
  thetamin = theta(Mminindex,:);
  
  for i=1:nD,
    thetamine(i,:) = (theta(:,i) - (1/Nmin)*sum(Omin,2))'*Zc*bay.Rscores;
  end

  term1 = kxx - 2/(Nmin)*sum(thetamin,1)';
  term2 = (Nmin^-2)*sum(sum(Ominmin));
  term3 = thetamine.^2*Diagmatrix;
  
  var_min = (term1+term2)./bay.mu-term3;
  
  
% Ppos, Pmin, res
for i=1:nD,    
  pdfplus = priorpos   * normpdf(Py(i),Pymp,sqrt(1/bay.zeta+var_plus(i)));
  pdfmin = (1-priorpos)* normpdf(Py(i),Pymn,sqrt(1/bay.zeta+var_min(i)));
  som = pdfmin+pdfplus;
  Pplus(i,1) = pdfplus./som;
  Pmin(i,1) = pdfmin./som;
end