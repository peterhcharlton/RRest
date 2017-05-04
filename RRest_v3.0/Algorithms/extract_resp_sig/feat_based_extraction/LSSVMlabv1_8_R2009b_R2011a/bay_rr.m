function [A,B,C,D,E,F,G] = bay_rr(X,Y,gam,level,nb,eigvals,eigvec)
% Bayesian inference of the cost on the three levels of linear ridge regression
% 
% >> cost = bay_rr(X, Y, gam, level)
% 
% This function implements the cost functions related to the
% Bayesian framework of linear ridge Regression [29]. Optimizing
% this criteria results in optimal model parameters W,b,
% hyperparameters. The criterion can also be used for model
% comparison. 
% 
% The obtained model parameters w and b are optimal on the first
% level w.r.t J = 0.5*w'*w+gam*0.5*sum(Y-X*w-b).^2. 
% 
% Full syntax
% 
%     * OUTPUTS on the first level: Cost proportional to the posterior of the model parameters.
% 
% >> [costL1, Ed, Ew] = bay_rr(X, Y, gam, 1)
%   
%         costL1: Cost proportional to the posterior
%         Ed(*) : Cost of the fitting error term
%         Ew(*) : Cost of the regularization parameter
% 
%     * OUTPUTS on the second level: Cost proportional to the posterior of gam.
% 
% >> [costL2, DcostL2, Deff, mu, ksi, eigval, eigvec] = bay_rr(X, Y, gam, 2)
% 
%         costL2    : Cost proportional to the posterior on the second level
%         DcostL2(*): Derivative of the cost proportional to the posterior
%         Deff(*)   : Effective number of parameters
%         mu(*)     : Relative importance of the fitting error term
%         ksi(*)    : Relative importance of the regularization parameter
%         eigval(*) : Eigenvalues of the covariance matrix
%         eigvec(*) : Eigenvectors of the covariance matrix
% 
%     * OUTPUTS on the third level: The following commands can be
%     used to compute the level 3 cost function for different
%     models (e.g. models with different selected sets of
%     inputs). The best model can then be chosen as the model with
%     best level 3 cost (CostL3). 
% 
% >> [costL3, gam_optimal] = bay_rr(X, Y, gam, 3)
% 
%         costL3         : Cost proportional to the posterior on the third inference level
%         gam_optimal(*) : Optimal regularization parameter obtained from optimizing the second level
% 
%     * INPUTS:
% 
% >> cost = bay_rr(X, Y, gam, level)
% 
%         X N x d matrix with the inputs of the training data
%         Y N x 1 vector with the outputs of the training data
%         gam Regularization parameter
%         level 1, 2, 3
% 
% See also:
%   ridgeregress,bay_lssvm

% -dislaimer-

  
  if ~exist('fminunc'),
    error('This function needs the optimization function ''fminunc''.');
  end    
  
[N,d] = size(X);

if (level==1 || level==2),
  [W,b] = ridgeregress(X,Y,gam);
  Ew =  .5*W'*W;             % Ew
  Ed = .5*sum((X*W+b-Y).^2); % Ed
  Ewgd = Ew+gam*Ed;          % Ewgd
  A=Ewgd; C=Ed; B=Ew;
  
  if level==2,
    if nargin>=7,
      if numel(eigvals)>length(eigvals), v = diag(eigvals); else v=eigvals; end
      V = eigvec;      
    else
      eval('[V,v] = eigs(X''*X+eye(d)*2,nb);','[V,v] = eig(X''*X+eye(d)*2);v=diag(v);');v=v-2;
      v = v*(N-1)/N;
    end
    Peff = find(v>1000*eps); Neff=length(Peff);
    v = v(Peff); V= V(:,Peff); 
    vall = zeros(N-1,1);vall(1:Neff)=v;
    Deff = 1+sum(gam.*v./(1+gam.*v));                 % Deff
    CostL2 = sum(log(vall+1./gam))+(N-1)*log(Ewgd);   % CostL2
    DcostL2 = -sum(1./(gam+vall.*gam^2))+(N-1)*Ed/Ewgd; % DcostL2
    mu = 2*Ed/(N-Deff); ksi = mu*gam;                 % mu and ksi
    A=CostL2; B=DcostL2; C=Deff; D=mu; E=ksi;F=v;G=V;
  end 
 
elseif level==3,
  
  % check fminunc
  resp = which('fminunc');
  %disp(' ');
  if isempty(resp),
    error(' ''fminunc'' not available');
  end

  eval('nb;','nb=''blabla'';');
  opties=optimset('MaxFunEvals', 2000,'GradObj','on', 'DerivativeCheck', 'off', 'TolFun', .0001, 'TolX', .0001, 'Display','off' );
  [CostL2,DCostL2, Deff, mu,ksi,v, V] = bay_rr(X,Y,gam,2,nb);
  gam_opt = exp(fminunc(@costL2, log(gam), opties,X,Y,nb,v,V));
  [CostL2,DCostL2, Deff, mu,ksi,v, V] = bay_rr(X,Y,gam_opt,2,nb,v,V);  
  CostL3 = .5*length(v)*log(mu)+(N-1)*log(ksi) - log(Deff-1)-log(N-Deff) - sum(log(mu+ksi*v));
  A = CostL3; B = gam_opt;

else
  error('level should be ''1'', ''2'' or ''3''.');
end




function [C,Dc] = costL2(log_gam,X,Y,nb,v,V)
%
[C,Dc] = bay_rr(X,Y,exp(log_gam),2,nb,v,V); 
