function [cost,retained] = trimmedmse(R,beta,V);
% Calculate trimmed mean of the squared value of the residuals.
%
%  cost = trimmedmse(R);
%
% The factor where one trimms off the normed residuals is
% optimized. However, one can pass a default value when one to
% exclude this optimization, e.g.:
%
%  cost = trimmedmse(R,0.15);
%
% One can overrule the default norm (norm='abs') by passing the norm function.
% 
%  cost = trimmedmse(R,[],norm);
%
% As an additional output, the index of the retained points ca be
% received:
%
%  [cost,retained] = trimmedmse(R);
%
% see also:
%   mse, misclass

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab



% default trimming?
eval('beta;','beta=[];');
eval('R = feval(V,R);','R = R.^2;');
[Rs,si] = sort(R);
N = max(size(Rs));

%figure; hist(Rs,50); pause

if ~isempty(beta),
  nb = N - floor(N*beta);
  mu = mean(Rs(1:nb));  
  cost = mu;
else
% optimize trimming factor  
  best_variance = inf;
  t = 1;
  %betas = 0:.01:.45;
  %betas = [0 0.05 0.10 0.175 0.30 0.45];
  betas = 0.05;
  for beta = betas,
    
    %beta = beta*2;
    nb = N - floor(N*beta);
    mu = mean(Rs(1:nb));  
    %variance = 1/((1-beta)^2) * (sum((Rs(1:nb)-mu).^2)/N+ (beta*(Rs(nb)-mu)^2));
    variance = sum((Rs(1:nb)-mu).^2) + ...
        (floor(N*beta)+1)*(Rs(nb)-mu)^2 - ...
        1/N*(floor(N*beta)*(Rs(nb)-mu)^2);
    variance = variance/(nb*nb-1);
    %v(t,1) = variance; t=t+1;
    if variance <= best_variance,
      best_variance = variance;
      cost = mu;
      best_beta = beta;
    end
  end
end
%figure; plot(betas',[v sum(v,2)]);
%figure; hist(Rs,50); 
% which are the retained data points
retained = si(1:nb);



