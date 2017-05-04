% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


clc;

disp('This is a simple demo, solving a simple regression task using');
disp('LS-SVMlab. A dataset is constructed in the right formatting. The');
disp('data are represented as matrices where each row contains one');
disp('datapoint: ');
disp(' ');
disp('press <ENTER> key'); pause
disp(' ');

disp('>> X = (-3:0.2:3)'';');
X = (-3:0.2:3)';
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
disp('parameter of the RBF kernel:');

disp(' ');
disp('>> gam = 10;');
gam = 10;
disp('>> sig2 = 0.3;');
sig2 = 0.3;
disp('>> type = ''function estimation'';');
type = 'function estimation';
disp('>> [alpha,b] = trainlssvm({X,Y,type,gam,sig2,''RBF_kernel''});');
[alpha,b] = trainlssvm({X,Y,type,gam,sig2,'RBF_kernel'});
disp(' ');
disp('press <ENTER> key'); pause
disp(' ');

disp('The parameters and the variables relevant for the LS-SVM are');
disp('passed as one cell. This cell allows for consistent default');
disp('handling of LS-SVM parameters and syntactical grouping of related');
disp('arguments. This definition should be used consistently throughout');
disp('the use of that specific LS-SVM model');
disp('The object oriented interface to LS-SVMlab leads to');
disp('shorter function calls (see demomodel). ');

disp('By default, the data are preprocessed by application of the function');
disp('prelssvm to the raw data and the function postlssvm on the');
disp('predictions of the model. This option can explicitly be switched off in');
disp('the call: ');

disp(' ');
disp('>> [alpha,b] = trainlssvm({X,Y,type,gam,sig2,''RBF_kernel'',''original''});');
[alpha,b] = trainlssvm({X,Y,type,gam,sig2,'RBF_kernel','original'});
disp(' ');

disp('or can be switched on (default):');

disp(' ');
disp('>> [alpha,b] = trainlssvm({X,Y,type,gam,sig2,''RBF_kernel'',''preprocess''});');
[alpha,b] = trainlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'});

%disp('Remember to consistently use the same option in all successive calls.');
disp(' ');
disp('press <ENTER> key'); pause
disp(' ');

disp('To evaluate new points for this model, the function simlssvm is');
disp('used. At first, test data is generated: ');

disp(' ');
disp('>> Xt = 3.*randn(10,1);');
Xt = 3.*randn(10,1);
disp(' ');

disp('Then, the obtained model is simulated on the test data:'); 

disp(' ');
disp('>> Yt = simlssvm({X,Y,type,gam,sig2,''RBF_kernel'',''preprocess''},{alpha,b},Xt);');
Yt = simlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b},Xt);
disp(' ');
disp('>> Y');
Y
disp(' ');
disp('press <ENTER> key'); pause
disp(' ');


disp('The LS-SVM result can be displayed if the dimension of the input');
disp('data is 1 or 2. ');

disp(' ');
disp('>> plotlssvm({X,Y,type,gam,sig2,''RBF_kernel'',''preprocess''},{alpha,b});');
figure; plotlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b});
disp(' ');

disp('All plotting is done with this simple command. It looks for the');
disp('best way of displaying the result. When the real function is known,');
disp('it can be displayed as follows:');
disp(' ');
disp('>> hold on; plot(min(X):.1:max(X),sinc(min(X):.1:max(X)),''r-.'');'); hold off
Xt = (min(X):.1:max(X))'; 
eval('Yt = sinc(Xt);',...
     'Yt = sin(pi.*Xt+12345*eps)./(pi*Xt+12345*eps)+0.1.*randn(length(Xt),1);');
hold on;  plot(Xt,Yt,'r-.'); hold off
disp(' ');
disp(' This concludes the demo');