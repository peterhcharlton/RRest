% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


clc;

disp('This is a simple demo, solving a simple regression task using');
disp('LS-SVMlab and constructing confidence intervals. A dataset is constructed in the right formatting. The');
disp('data are represented as matrices where each row contains one');
disp('datapoint: ');
disp(' ');
disp('press <ENTER> key'); pause
disp(' ');

disp('>> X = (-3:0.02:3)'';');
X = (-3:0.02:3)';
disp('>> Y = sinc(X)+0.1.*randn(length(X),1);');
eval('Y = sinc(X)+0.1.*randn(length(X),1);',...
     'Y = sin(pi.*X+12345*eps)./(pi*X+12345*eps)+0.1.*randn(length(X),1);');
disp('>> X');
X

disp('>> Y');
Y

disp('In order to make an LS-SVM model, we need 2 extra parameters: gamma');
disp('(gam) is the regularization parameter, determining the trade-off');
disp('between the fitting error minimization and smoothness of the');
disp('estimated function. sigma^2 (sig2) is the kernel function');
disp('parameter of the RBF kernel. These can be found via cross-validation:');
disp(' ');

model = initlssvm(X,Y,'f',[],[],'RBF_kernel','o');

disp('>> model = tunelssvm(model,''simplex'',''crossvalidatelssvm'',{10,''mse''});');

model = tunelssvm(model,'simplex','crossvalidatelssvm',{10,'mse'});
disp(' ');
disp('press <ENTER> key'); pause
disp(' ');

disp('Training the model ');

disp(' ');
disp('>> model = trainlssvm(model)');
model = trainlssvm(model);
disp(' ');



disp('Computation of Confidence Intervals ');
disp(' ');
disp('press <ENTER> key'); pause


disp('ci = cilssvm(model);')
ci = cilssvm(model);


disp('The LS-SVM result and confidence intervals can be displayed if the dimension of the input');
disp('data is 1 or 2. ');

disp(' ');
disp('>> plotlssvm(model);');
figure; plotlssvm(model);
disp(' ');

hold all
fill([X;flipud(X)],[ci(:,1);flipud(ci(:,2))],'c','FaceAlpha',0.5,'EdgeAlpha',1,'EdgeColor','w')

disp('All plotting is done with this simple command. It looks for the');
disp('best way of displaying the result.') 
disp(' ');
disp(' This concludes the demo');