% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


clc;

disp('A simple example shows how to start using the toolbox for a');
disp('classification task. We start with constructing a simple example');
disp('dataset according to the right formatting. Data are represented ');
disp('as matrices where each row contains one datapoint: ');
disp(' ');
disp('press <ENTER> key'); pause
disp(' ');

disp(' >> X = 2.*rand(30,2)-1;');
X = 2.*rand(30,2)-1;
disp(' >> Y = sign(sin(X(:,1))+X(:,2));');
Y = sign(sin(X(:,1))+X(:,2));
disp(' >> X');
X

disp(' >> Y');
Y

disp(' ');
disp('press <ENTER> key'); pause
disp(' ');

disp('In order to make an LS-SVM model, we need 2 extra parameters: gamma');
disp('(gam) is the regularization parameter, determining the trade-off');
disp('between the fitting error minimization and smoothness. In the');
disp('common case of the RBF kernel, sigma^2 (sig2) is the bandwidth:');
disp(' ');
disp(' >> gam = 10;');
gam = 10;
disp(' >> sig2 = 0.2;');
sig2 = 0.2;
disp(' >> type = ''classification'';');
type = 'classification';
disp(' >> [alpha,b] = trainlssvm({X,Y,type,gam,sig2,''RBF_kernel''});');
[alpha,b] = trainlssvm({X,Y,type,gam,sig2,'RBF_kernel'});

disp(' ');
disp('press <ENTER> key'); pause
disp(' ');

disp('The parameters and the variables relevant for the LS-SVM are');
disp('passed as one cell. This cell allows for consistent default');
disp('handling of LS-SVM parameters and syntactical grouping of related');
disp('arguments. This definition should be used consistently throughout');
disp('the use of that specific LS-SVM model.');
disp('The corresponding object oriented interface');
disp('to LS-SVMlab leads to shorter function calls (see demomodel). ');

disp('By default, the data are preprocessed by application of the function');
disp('prelssvm to the raw data and the function postlssvm on the');
disp('predictions of the model. This option can explicitly be switched off in');
disp('the call: ');

disp(' ');
disp(' >> [alpha,b] = trainlssvm({X,Y,type,gam,sig2,''RBF_kernel'',''original''});');
[alpha,b] = trainlssvm({X,Y,type,gam,sig2,'RBF_kernel','original'});
disp(' ');
disp('or be switched on (by default):');
disp(' ');
disp(' >> [alpha,b] = trainlssvm({X,Y,type,gam,sig2,''RBF_kernel'',''preprocess''});');
[alpha,b] = trainlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'});

disp(' ');
disp('press <ENTER> key'); pause
disp(' ');

%disp('Remember to consistently use the same option in all successive calls');
disp('To evaluate new points for this model, the function');
disp('simlssvm is used:');
disp(' ');
disp(' >> Xt = 2.*rand(10,2)-1;');
Xt = 2.*rand(10,2)-1;
disp(' >> Ytest = simlssvm({X,Y,type,gam,sig2,''RBF_kernel'',''preprocess''},{alpha,b},Xt);');
Ytest = simlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b},Xt);
disp(' ');
disp('The LS-SVM result can be displayed if the dimension of the input');
disp('data is 2. ');

disp(' ');
disp('press <ENTER> key'); pause
disp(' ');

disp(' ');
disp(' >> plotlssvm({X,Y,type,gam,sig2,''RBF_kernel'',''preprocess''},{alpha,b});');
figure; plotlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b});
disp(' ');
disp('All plotting is done with this simple command. It looks for the');
disp('best way of displaying the result. ');
disp(' ');
disp(' This concludes the demo');
