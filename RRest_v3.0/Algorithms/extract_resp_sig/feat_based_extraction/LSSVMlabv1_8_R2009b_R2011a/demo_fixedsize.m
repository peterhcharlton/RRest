% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


disp(' This demo illustrates the idea of fixed size LS-SVM. ');
disp(' The program consists of 2 steps. In the former, one ');
disp(' constructs a reduced set of support vectors base on ');
disp(' an apropriate criterion on the data. In this case');
disp(' the measure ''kentropy'' is optimized.');
disp(' ');
disp(' In the latter step, one constructs the implicit mapping');
disp(' to feature space based on the eigenvalue decomposition.');
disp(' A parametric linear regression is executed on the mapped');
disp(' data');
disp(' ');
disp(' To see the used cose, use the call');
disp(' ');
disp('>> type demo_fixedsize ');
disp(' ');
disp(' or ');
disp(' ');
disp('>> edit demo_fixedsize ');
disp(' ');
disp(' A dataset is constructed at first...');
%
% dataset
%clear
figure;
randn('state',0);
x = sort(2.*randn(2000,1));
x0 = sort(2.*randn(2000,1));
%x=(-9.95:0.1:10)';
%x0=(-9.95:0.05:10)';
eval('y = sinc(x)+0.05.*randn(length(x),1);',...
     'y = sin(pi.*x+12345*eps)./(pi*x+12345*eps)+0.05.*randn(length(x),1);');

eval('y0 = sinc(x0)+0.05.*randn(length(x0),1);',...
     'y0 = sin(pi.*x0+12345*eps)./(pi*x0+12345*eps)+0.05.*randn(length(x0),1);');

disp(' The parameters are initialized...');


%
% initiate values
kernel = 'RBF_kernel';
sigma2=.75;
gamma=1;
crit_old=-inf;
Nc=15;
Xs=x(1:Nc,:);
Ys=y(1:Nc,:);

disp(' The optimal reduced set is constructed iteratively: ');

%
% iterate over data
%
tv = 1;
for tel=1:length(x)
   
  
  %
  % new candidate set
  %
  Xsp=Xs; Ysp=Ys;
  S=ceil(length(x)*rand(1));
  Sc=ceil(Nc*rand(1));
  Xs(Sc,:) = x(S,:);
  Ys(Sc,:) = y(S);
  Ncc=Nc;

  %
  % automaticly extract features and compute entropy
  %
  crit = kentropy(Xs,kernel, sigma2);
  
  if crit <= crit_old,
    crit = crit_old;
    Xs=Xsp;
    Ys=Ysp;
  else
    crit_old = crit;

    %
    % ridge regression    
    %
    [features,U,lam] = AFEm(Xs,kernel, sigma2,x);
    [w,b,Yh] = ridgeregress(features,y,gamma,features);

    %
    % make-a-plot
    %
    plot(x,y,'*'); hold on
    plot(x,Yh,'r-')
    plot(Xs,Ys,'go','Linewidth',7)
    xlabel('X'); ylabel('Y'); 
    title(['Approximation by fixed size LS-SVM based on maximal entropy: ' num2str(crit)]);
    hold off;  drawnow
  
   
  end
  
  %
  % validate
  %
  %Yh0 = AFE(Xs,kernel, sigma2,x0)*w + b;
  %val(tv,2) = mse(Yh0-y0); tv=tv+1;
  
end



disp(' The parametric linear ridge regression is calculated:');

%
% ridge regression    
%
features = AFEm(Xs,kernel, sigma2,x);    

% Bayesian inference of the gamma
try,
  [CostL3, gamma_optimal] = bay_rr(features,y,gamma,3);
catch,
  warning('no Bayesian optimization of the regularization parameter');
  gamma_optimal = gamma;
end

[w,b] = ridgeregress(features,y,gamma_optimal);
Yh0 = AFEm(Xs,kernel, sigma2,x0)*w+b;
echo off;         

%
% make-a-plot
plot(x,y,'*'); hold on
plot(x0,Yh0,'r-')
plot(Xs,Ys,'go','Linewidth',7)
xlabel('X'); ylabel('Y'); 
title(['Approximation by fixed size LS-SVM based on maximal entropy: ' num2str(crit)]);
hold off;  





