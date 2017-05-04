% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

%Multi-class example
clear all 
clc
disp('This is a simple demo, solving a simple multiclass classification problem with');
disp('LS-SVMlab using the object interface. The problem is solved using the One vs. One coding scheme.');
disp('A dataset is constructed in the right formatting.') 
disp(' ');
disp('press <ENTER> key'); pause
disp(' ');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate multi-class data composed out of a mixture of 2 gaussians

DIM=2;SIZE=50;SIZEte=5000;SIZE=floor(SIZE/2)*2;SIZEte=floor(SIZEte/2)*2;

randn('state',0);
X=[];
X=[X ; 1.6*randn(SIZE/2,DIM)+repmat([0 0],SIZE/2,1) ];
X=[X ; 0.9*randn(SIZE/2,DIM)+repmat([2 0],SIZE/2,1) ]; 
X=[X ; 0.8*randn(SIZE/2,DIM)+repmat([-1 1],SIZE/2,1)];
X=[X ; 0.9*randn(SIZE/2,DIM)+repmat([-1.3 3.5],SIZE/2,1)];
X=[X ; 1*randn(SIZE/2,DIM)+repmat([-2 1],SIZE/2,1)];
X=[X ; 0.9*randn(SIZE/2,DIM)+repmat([-3.5 0.2],SIZE/2,1)] ;
%[X,m,v,s]=standardize(X);
y=[];for i=1:3, y= [y ; i*ones(SIZE,1)]; end

Xt=[];
Xt=[Xt ; 1.6*randn(SIZEte/2,DIM)+repmat([0 0],SIZEte/2,1) ];
Xt=[Xt ; 0.9*randn(SIZEte/2,DIM)+repmat([2 0],SIZEte/2,1) ]; 
Xt=[Xt ; 0.8*randn(SIZEte/2,DIM)+repmat([-1 1],SIZEte/2,1)];
Xt=[Xt ; 0.9*randn(SIZEte/2,DIM)+repmat([-1.3 3.5],SIZEte/2,1)];
Xt=[Xt ; 1*randn(SIZEte/2,DIM)+repmat([-2 1],SIZEte/2,1)];
Xt=[Xt ; 0.9*randn(SIZEte/2,DIM)+repmat([-3.5 0.2],SIZEte/2,1)] ;

%Xt=standardize(Xt,m,v,s);
yt=[];for i=1:3, yt= [yt ; i*ones(SIZEte,1)];end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
disp('Start initialization, tuning procedure and training');
disp('');
disp('>> model = initlssvm(X,y,''c'',[],[],''RBF_kernel'');')
disp('>> model = tunelssvm(model,''simplex'',''crossvalidatelssvm'',{10,''misclass''},''code_OneVsOne'');')
disp('>> model = trainlssvm(model);')
disp('press <ENTER> key'); pause
disp(' ');

t1=cputime;
model = initlssvm(X,y,'c',[],[],'RBF_kernel');
model = tunelssvm(model,'simplex','crossvalidatelssvm',{10,'misclass'},'code_OneVsOne');
model = trainlssvm(model);
Y = simlssvm(model,Xt);

t2=cputime;
fprintf(1,'Tuning time %i \n',t2-t1);
fprintf(1,'Accuracy: %2.2f\n',100*sum(Y==yt)/length(yt));
plotlssvm(model,[],150);


