function m = tbform(model,alpha)

%
% INTERNAL FUNCTION
%

% This function calculates the factor 'm' for simultaneous CI/PI using the
% tube formula based "m=sqrt(2log(kappa/(alpha*pi)))" on upcrossing theory:
%
% Rice, S. O. (1939). The distribution of the maxima of a random curve.
% American Journal of Mathematics 61: 409–16.

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


if nargin <= 1
    alpha = 0.05;
end

% version check
resp = which('quadgk');
if isempty(resp), imethod=@quad; else imethod=@quadgk; end

% set options for 'fsolve' and 'fzero'
opt = optimset('Display','off','TolX',1e-15);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some initial calculation used later
K = kernel_matrix(model.xtrain,model.kernel_type,model.kernel_pars);
%Z = pinv(K+eye(model.nb_data)./model.gam);
Z = pinv(K+diag(ones(model.nb_data,1)./model.gam(:)));
c = sum(sum(Z));
J = ones(model.nb_data)./c;
S = K*(Z-Z*J*Z) + J*Z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(model.xtrain,2) == 1 && model.preprocess(1)=='o'
    
    % Calculation of the kappa coefficient in the tube formula
    kappa = imethod(@(x)fun1(model,x,Z,J),min(model.xtrain),max(model.xtrain));
    
    % Tube formula
    df = model.nb_data - trace(S);
    tube = @(x)(2*(1-tcdf(x,df))+(kappa/pi)*(1+(x.^2/df)).^(-df/2))-alpha;
    m = fsolve(tube,3,opt);    
    
else % higher dimensions (QMC integration)
    latticeseq_b2 init0
    e = 256;
    O = latticeseq_b2(size(model.xtrain,2), e)';
    O = latticeseq_b2(size(model.xtrain,2), e)';
    [O,dom] = scale(O,model.xtrain);
    for i=1:e, I(:,1) = kfunmore(model,O(i,:),Z,J);end
    kappa = mean(I)*dom;
    xi0 = (pi/2)*(trace(S(1:end-1,1:end-1))-1);
    d = size(model.xtrain,2);
    W = (eye(model.nb_data)-S')*(eye(model.nb_data)-S);
    df = (trace(W))^2/trace(W^2);
    tube = @(x)kappa*(gamma((d+1)/2)./pi^((d+1)/2)).*(1-fcdf(x.^2/(d+1),d+1,df))+ (xi0/2).*(gamma(d/2)./pi^(d/2)).*(1-fcdf(x.^2/d,d,df))-alpha;
    m = fzero(tube,3,opt);
    %m = sqrt(2*log(kappa/(pi*alpha)));
end


function I = fun1(model,x,Z,J)
L = smootherlssvm(model,x')';
dL = fund(model,x',Z,J)';
I = sqrt(sum(L.^2).*sum(dL.^2) - diag(L'*dL)'.^2)./sum(L.^2);%diag(sqrt((h.*(dL*dL'))-(L*dL').^2)./h)';

function I = kfunmore(model,x,Z,J)
L = smootherlssvm(model,x)';
dL = fund(model,x,Z,J);
[~,R] = qr([L dL],0);
I = prod(diag(R))/R(1,1)^2;

function dL = fund(model,x,Z,J)
% Calculation of the elementwise derivites (analitically) of the
% smoother matrix
Kt = kernel_matrix(model.xtrain,model.kernel_type,model.kernel_pars,x);
if size(x,2) == 1
    omega = repmat(x,1,model.nb_data);
    omega = omega - repmat(model.xtrain',size(omega,1),1);
    dL = -((1/model.kernel_pars)*omega.*Kt')*(Z-Z*J*Z);
else
    omega = repmat(x,size(model.xtrain,1),1)-model.xtrain;
    dL = -(1/model.kernel_pars)*(Z-Z*J*Z)'*(omega.*repmat(Kt,1,size(x,2)));
end

function [O,dom] = scale(O,x)
mi = min(x);
ma = max(x);
dom = prod(ma-mi);
for i=1:length(mi)
    O(:,i) = mi(i) + (ma(i)-mi(i))*O(:,i);
end