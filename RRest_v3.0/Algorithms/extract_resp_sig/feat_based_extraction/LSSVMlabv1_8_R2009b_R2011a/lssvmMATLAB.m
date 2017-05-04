function [model,H] = lssvmMATLAB(model) 
% Only for intern LS-SVMlab use;
%
% MATLAB implementation of the LS-SVM algorithm. This is slower
% than the C-mex implementation, but it is more reliable and flexible;
%
%
% This implementation is quite straightforward, based on MATLAB's
% backslash matrix division (or PCG if available) and total kernel
% matrix construction. It has some extensions towards advanced
% techniques, especially applicable on small datasets (weighed
% LS-SVM, gamma-per-datapoint)

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


%fprintf('~');
%
% is it weighted LS-SVM ?
%
weighted = (length(model.gam)>model.y_dim);
if and(weighted,length(model.gam)~=model.nb_data),
  warning('not enough gamma''s for Weighted LS-SVMs, simple LS-SVM applied');
  weighted=0;
end

% computation omega and H
omega = kernel_matrix(model.xtrain(model.selector, 1:model.x_dim), ...
    model.kernel_type, model.kernel_pars);


% initiate alpha and b
model.b = zeros(1,model.y_dim);
model.alpha = zeros(model.nb_data,model.y_dim);

for i=1:model.y_dim,
    H = omega;
    model.selector=~isnan(model.ytrain(:,i));
    nb_data=sum(model.selector);
    if size(model.gam,2)==model.nb_data, 
      try invgam = model.gam(i,:).^-1; catch, invgam = model.gam(1,:).^-1;end
      for t=1:model.nb_data, H(t,t) = H(t,t)+invgam(t); end
    else
      try invgam = model.gam(i,1).^-1; catch, invgam = model.gam(1,1).^-1;end
      for t=1:model.nb_data, H(t,t) = H(t,t)+invgam; end
    end    

    v = H(model.selector,model.selector)\model.ytrain(model.selector,i);
    %eval('v  = pcg(H,model.ytrain(model.selector,i), 100*eps,model.nb_data);','v = H\model.ytrain(model.selector, i);');
    nu = H(model.selector,model.selector)\ones(nb_data,1);
    %eval('nu = pcg(H,ones(model.nb_data,i), 100*eps,model.nb_data);','nu = H\ones(model.nb_data,i);');
    s = ones(1,nb_data)*nu(:,1);
    model.b(i) = (nu(:,1)'*model.ytrain(model.selector,i))./s;
    model.alpha(model.selector,i) = v(:,1)-(nu(:,1)*model.b(i));
end
return




