% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

disp('  This demo illustrates facilities of LS-SVMlab');
disp('  with respect to unsupervised learning.');

disp(' a demo dataset is generated...');
clear yin yang samplesyin samplesyang mema
% initiate variables and construct the data
nb =200;
sig = .20;

% construct data
leng = 1;
for t=1:nb, 
  yin(t,:) = [2.*sin(t/nb*pi*leng) 2.*cos(.61*t/nb*pi*leng) (t/nb*sig)]; 
  yang(t,:) = [-2.*sin(t/nb*pi*leng) .45-2.*cos(.61*t/nb*pi*leng) (t/nb*sig)]; 
  samplesyin(t,:)  = [yin(t,1)+yin(t,3).*randn   yin(t,2)+yin(t,3).*randn];
  samplesyang(t,:) = [yang(t,1)+yang(t,3).*randn   yang(t,2)+yang(t,3).*randn];
end

% plot the data
figure; hold on;
plot(samplesyin(:,1),samplesyin(:,2),'+','Color',[0.6 0.6 0.6]);
plot(samplesyang(:,1),samplesyang(:,2),'+','Color',[0.6 0.6 0.6]);
xlabel('X_1');
ylabel('X_2');
title('Structured dataset');
disp('  (press any key)');
pause

%
% kernel based Principal Component Analysis
%
disp(' ');
disp('  extract the principal eigenvectors in feature space');
disp(' >> nb_pcs=4;'); nb_pcs = 4;
disp(' >> sig2 = .8;'); sig2 = .8;
disp(' >> [lam,U] = kpca([samplesyin;samplesyang],''RBF_kernel'',sig2,[],''eigs'',nb_pcs); ');
[lam,U] = kpca([samplesyin;samplesyang],'RBF_kernel',sig2,[],'eigs',nb_pcs);
disp('  (press any key)');
pause

%
% make a grid over the inputspace
%
disp(' ');
disp(' make a grid over the inputspace:');
disp('>> Xax = -3:0.1:3; Yax = -2.0:0.1:2.5;'); Xax = -3:0.1:3; Yax = -2.0:0.1:2.5;
disp('>> [A,B] = meshgrid(Xax,Yax);'); [A,B] = meshgrid(Xax,Yax);
disp('>> grid = [reshape(A,prod(size(A)),1) reshape(B,1,prod(size(B)))'']; ');
grid = [reshape(A,numel(A),1) reshape(B,1,numel(B))'];


%
% compute projections of each point of the inputspace on the
% principal components
%
disp(' ');
disp('  compute projections of each point of the inputspace on the  ');
disp('  principal components');
disp('>> k = kernel_matrix([samplesyin;samplesyang],''RBF_kernel'',sig2,grid)'';  ');
k = kernel_matrix([samplesyin;samplesyang],'RBF_kernel',sig2,grid)';
disp('>> projections = k*U;'); projections = k*U; 
disp('>> contour(Xax,Yax,reshape(projections(:,1),length(Yax),length(Xax)));'); 
contour(Xax,Yax,reshape(projections(:,1),length(Yax),length(Xax)));
title('Projections onto the first kernel PC');
disp('  (press any key)');
pause



%
% Compute the approximate pre-image in the input space
disp(' ');
disp(' Compute the approximate pre-image in the input space');


  disp(' For every point, the approximate pre-image is computed using:');
  disp(' ----------------------------------------------------------');
  disp(' ');
  disp('>> Xd=preimage_rbf([samplesyin;samplesyang],sig2,U); ');
  Xd=preimage_rbf([samplesyin;samplesyang],sig2,U);
  figure; hold on;
  plot(samplesyin(:,1),samplesyin(:,2),'+','Color',[0.6 0.6 0.6]);
  plot(samplesyang(:,1),samplesyang(:,2),'+','Color',[0.6 0.6 0.6]);
  xlabel('x_1');
  ylabel('x_2');
  disp('>> plot(Xd(:,1),Xd(:,2),''ko''); '); plot(Xd(:,1),Xd(:,2),'bo');
  disp(' ');
  title('Denoising (''o'') by computing an approximate pre-image');
  disp(' ');
  disp(' In the last figure, one can see the original datapoints');
  disp('(''*'') and the reconstructed data (''o'').  ');
  disp(' ');
  disp(' ');
disp('  This  concludes this demo');
hold off 