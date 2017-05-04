% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

disp(' This demo illustrates how one can use different components ');
disp(' of the toolbox to compose a fixed size LS-SVM classifier.');
disp(' To study how this result can be achieved, study the used commands.');
disp(' These can be shown by the command');
disp(' ');
disp(' >> type demo_fixedclass');
disp(' ');
disp(' or ');
disp(' ');
disp(' >> edit demo_fixedclass');
disp(' ');


load ripley

%
% initiate values
type = 'classification';
gamma = 0.1;
kernel = 'RBF_kernel';
sigma2 = 1;
sigma2ent = 0.1;
crit_old=-inf;
Nc=20;
Xs=X(1:Nc,:);
Ys=Y(1:Nc,:);

%
% Initiate grid for plot
grain = 25;
xmin1=min(X(:,1)); 
xmax1=max(X(:,1)); 
xmin2=min(X(:,2)); 
xmax2=max(X(:,2)); 
xrange1 = xmin1:(xmax1-xmin1)/grain:xmax1;
xrange2 = xmin2:(xmax2-xmin2)/grain:xmax2;
[XX,YY] = meshgrid(xrange1,xrange2);
Xt = [reshape(XX,numel(XX),1) reshape(YY,numel(YY),1)];
figure;


%
% iterate over data
%
for tel=1:5*length(X)
   
  
  %
  % new candidate set
  %
  Xsp=Xs; Ysp=Ys;
  S=ceil(length(X)*rand(1));
  Sc=ceil(Nc*rand(1));
  Xs(Sc,:) = X(S,:);
  Ys(Sc,:) = Y(S);
  Ncc=Nc;

  %
  % automaticly extract features and compute entropy
  %
  crit = kentropy(Xs,kernel, sigma2ent);
  
  if crit <= crit_old,
    crit = crit_old;
    Xs=Xsp;
    Ys=Ysp;
  else
    crit_old = crit;

    %
    % ridge regression    
    features   = AFEm(Xs,kernel, sigma2,X);
    features_t = AFEm(Xs,kernel, sigma2,Xt);
    [w,b,Yht] = ridgeregress(features,Y,gamma,features_t);
    Yht = sign(Yht);

    %
    % make-a-plot
    Ygt = reshape(Yht(:,1),size(XX,1),size(XX,2));
    colormap cool;
    [C,h]=contourf(XX,YY,Ygt); 
    hold on;
    n = find(Y<=0);
    np = plot(X(n,1),X(n,2),'k.'); 
    p = find(Y>0);
    pp = plot(X(p,1),X(p,2),'k+'); 
    sv = plot(Xs(:,1),Xs(:,2),'go','Linewidth',7);
    xlabel('X_1'); ylabel('X_2'); 
    title(['Approximation by fixed size LS-SVM based on maximal entropy: ' num2str(crit)]);
    legend([np pp sv],'Negative points','Positive points',...
           'Support Vectors');
    
    hold off;  drawnow
  
  end
    
end


