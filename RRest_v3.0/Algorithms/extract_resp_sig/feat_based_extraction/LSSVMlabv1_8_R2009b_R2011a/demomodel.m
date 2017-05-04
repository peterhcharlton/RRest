% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


clc;
disp(' This demo explains the use of the advanced object oriented interface');
disp('  ''model''. For first users, we recommend the functional interface');
disp('  as  explained in ''democlass'' and ''demofun''.');
disp('  ');
disp(' The ''model'' is the generic object which collects '); 
disp(' all relevant signals, parameters, options and functions ');
disp(' related to an application of the LS-SVM.');
disp(' This interface is depreciated for the casional users because');
disp(' of the implicit nature: the distinction between in- and ');
disp(' output vanishes. If one wants to use the full power  ');
disp(' of LS-SVMlab, one is recommended to go through this demo.');
disp(' ');
disp(' We focus on function estimation, however the insights are completely');
disp( 'equivalent for classification. A dataset is constructed at first:');
disp(' ');
disp('> X = (-3:.2:3)''');
X = (-3:.2:3)';
disp('> Y = X.^3+2.*randn(length(X),1);');
Y = X.^3+2.*randn(length(X),1);
disp(' ');
disp(' A model is CONSTRUCTED for this data');
disp(' ');
disp(' >> gam=1; sig2=1;'); gam=1; sig2=1;
disp(' >> model = initlssvm(X,Y,''function'',gam,sig2,''RBF_kernel'');');model = initlssvm(X,Y,'function',gam,sig2,'RBF_kernel');
disp(' ');
disp('press enter to continue...');
pause; 


disp('  The specifications of the model can be seen typing just the'); 
disp('  name of the object that is constructed.');
disp(' >> model');model
disp(' ');
disp(' If one wants to see the value of a specific option,');
disp(' the ''.'' operator is to be used:');
disp(' ');
disp(' >> model.preprocess'); model.preprocess
disp(' ');
disp('press enter to continue...');
pause; 

disp(' If one wants to CHANGE a value of a specific option of the');
disp(' model, the function ''changelssvm'' is to be used if consistent');
disp(' models are wanted:');
disp(' ');
disp(' >> model = changelssvm(model,''gam'',1.2);');model = changelssvm(model,'gam',1.2);
disp(' ');
disp(' The options can be divided in 4 classes, the general LS-SVM options,');
disp(' the trainpoint administration,');
disp(' the preprocess options and the encoding options.');
disp(' The help of ''changelssvm'' gives a description of the different');
disp(' fields of the model (type ''help changelssvm''). ');
disp(' The use of the model''s field ''status'' is important to see the full power');
disp(' of the model concept. The object ''model'' knows ');
disp(' whether it needs to be retrained. Retraining is needed if a specification');
disp(' of the model is changed since the last training.');
disp(' ');
disp(' The demo ... shows how to control the preprocessing and the coding');
disp(' using the appropriate model.');
disp(' ');
disp(' It is adviced to check carefully the model''s options before'); 
disp(' starting the calculations.');
disp(' ');
disp('press enter to continue...');
pause; 

disp(' If the model is clearly defined, the routine to train the model');
disp(' has to be called on this model. In the case of the LS-SVM, ');
disp(' this training is done by:');
disp(' ');
disp(' >> model = trainlssvm(model);');model = trainlssvm(model);
disp(' ');
disp(' Given the trained model, one can simulate some testpoints');
disp(' and make a plot of the model. (wait a few seconds...)');
disp(' ');
disp(' >> Xt = 2.*randn(10,1);');Xt = 2.*randn(10,1);
disp(' >> Yt = simlssvm(model,Xt)');Yt= simlssvm(model,Xt)
disp(' >> plotlssvm(model);');plotlssvm(model);
disp(' ');
disp(' By understanding this step, one masters basicly how to ');
disp(' use the object ''model''.');
disp(' ');
disp(' As an extra, the underlying function is also given as a dotted line')
disp(' on the plot');
disp('>> hold on; plot((-3:.2:3)'', sinc(-3:.2:3)'','':'');');
hold on; plot((-4:.2:4)', (-4:.2:4).^3','r:');
disp(' ');
disp('press enter to continue...');
pause; 


disp(' The same model can be used to try different options:');
disp(' ');
disp(' >> model = changelssvm(model,''gam'',3);'); model = changelssvm(model,'gam',3);
disp(' >> model = trainlssvm(model);'); model = trainlssvm(model);
disp(' >> plotlssvm(model);');plotlssvm(model);
disp(' ');
disp(' or with a polynomial kernel of degree 3');
disp(' ');
disp(' >> model = changelssvm(model,''kernel_type'',''poly_kernel'');'); model = changelssvm(model,'kernel_type','poly_kernel');
disp(' >> model = changelssvm(model,''kernel_pars'',[1;3]);');model = changelssvm(model,'kernel_pars',[1;3]);
disp(' >> model = trainlssvm(model);'); model = trainlssvm(model);
disp(' >> plotlssvm(model);');plotlssvm(model);
disp(' ');

disp(' This concludes this demo');

