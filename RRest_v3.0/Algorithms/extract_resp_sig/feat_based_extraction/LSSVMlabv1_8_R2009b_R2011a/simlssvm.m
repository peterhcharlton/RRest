function [Y,Yl,model] = simlssvm(model,Xt,A3,A4,A5)
% Evaluate the LS-SVM at the given points
%
% >> Yt = simlssvm({X,Y,type,gam,sig2,kernel}, {alpha,b}, Xt)
% >> Yt = simlssvm({X,Y,type,gam,sig2,kernel}, Xt)
% >> Yt = simlssvm(model, Xt)
%
% The matrix Xt represents the points one wants to predict. The
% first cell contains all arguments needed for defining the LS-SVM
% (see also trainlssvm, initlssvm). The second cell contains the
% results of training this LS-SVM model. The cell syntax allows for
% flexible and consistent default handling.
%
% As in training, three different implementations are included
% (simclssvm.mex*, simFILE.x, simlssvm.m). The cmex algorithm is
% called, except when specified otherwise. After a simulation call,
% a '.' is displayed.
%
%
% Full syntax
%
%     1. Using the functional interface:
%
% >> [Yt, Zt] = simlssvm({X,Y,type,gam,sig2}, Xt)
% >> [Yt, Zt] = simlssvm({X,Y,type,gam,sig2,kernel}, Xt)
% >> [Yt, Zt] = simlssvm({X,Y,type,gam,sig2,kernel,preprocess}, Xt)
% >> [Yt, Zt] = simlssvm({X,Y,type,gam,sig2,kernel,preprocess}, {alpha,b}, Xt)
%
%       Outputs
%         Yt            : Nt x m matrix with predicted output of test data
%         Zt(*)         : Nt x m matrix with predicted latent variables of a classifier
%       Inputs
%         X             : N x d matrix with the inputs of the training data
%         Y             : N x 1 vector with the outputs of the training data
%         type          : 'function estimation' ('f') or 'classifier' ('c')
%         gam           : Regularization parameter
%         sig2          : Kernel parameter (bandwidth in the case of the 'RBF_kernel')
%         kernel(*)     : Kernel type (by default 'RBF_kernel')
%         preprocess(*) : 'preprocess'(*) or 'original'
%         alpha(*)      : Support values obtained from training
%         b(*)          : Bias term obtained from training
%         Xt            : Nt x d inputs of the test data
%
%
%     2. Using the object oriented interface:
%
% >> [Yt, Zt, model] = simlssvm(model, Xt)
%
%       Outputs
%         Yt       : Nt x m matrix with predicted output of test data
%         Zt(*)    : Nt x m matrix with predicted latent variables of a classifier
%         model(*) : Object oriented representation of the LS-SVM model
%       Inputs
%         model    : Object oriented representation of the LS-SVM model
%         Xt       : Nt x d matrix with the inputs of the test data
%
% See also:
%   trainlssvm, initlssvm, plotlssvm, code, changelssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


%
% control inputs
%
if iscell(model),
    iscell_model = 1;
    model = initlssvm(model{:});
    if iscell(Xt),
        model.alpha = Xt{1};
        model.b = Xt{2};
        model.status = 'trained';
        eval('Xt = A3;',' ');
    end
    eval('nb_to_sim = A4;','nb_to_sim = size(Xt,1)-model.x_delays;');
    Yt = [];
else
    iscell_model = 0;
    if nargin>3,
        Yt = A3;
        eval('nb_to_sim = A4;','nb_to_sim = size(Xt,1)-model.x_delays;');
    else
        eval('nb_to_sim = A3;','nb_to_sim = size(Xt,1)-model.x_delays;');
        Yt =[];
    end
end

eval('Xt;','error(''Test data Xtest undefined...'');');

%
% check dimensions
%
if size(Xt,2)~=model.x_dim,
    error('dimensions of new datapoints Xt not equal to trainingsset...');
end
if ~isempty(Yt) && size(Yt,2)~=model.y_dim,
    error('dimensions of new targetpoints Yt not equal to trainingsset...');
end


%
% preprocessing ...
%
if model.preprocess(1)=='p',
    [Xt,Yt] = prelssvm(model,Xt,Yt);
end

%
% train if status is not  'trained'
%
if model.status(1)~='t', % not 'trained'
    warning('Model is not trained --> training now...')
    model = trainlssvm(model);
end

%
% if dimension of output >1
%
if model.y_dim>1,
    if length(model.kernel_type)>1 || size(model.kernel_pars,2)>1 || size(model.gam,2)==model.y_dim,
        %disp('multi dimensional output...');
        fprintf('m');
        [Y Yl] = simmultidimoutput(model,Xt,Yt,nb_to_sim);
        if iscell_model, model = Yl; end
        return
    end
end

%
% set parameters: how much points to evaluate and to simulate
%
if (model.type(1)=='c'),
    nb_sim = nb_to_sim;
    Yt=[];
elseif (model.type(1)=='f'),
    nb_sim = nb_to_sim;
    Yt=[];
end


%
% simulate the model (blockwize) using the MATLAB implementations;
%

bz=3000;
N=size(Xt,1);
NRofBlocks=floor(N/bz);
modu=N-NRofBlocks*bz;

Y=zeros(N,1);
for i=1:NRofBlocks,
    indb=(i-1)*bz+1:i*bz;
    Y(indb,:)=simFct(model,Xt(indb,:));
end

if modu~=0
    indb=NRofBlocks*bz+1:NRofBlocks*bz+modu;
    Y(indb,:)=simFct(model,Xt(indb,:));
end

%
% for classification
%
Yl = Y;
if model.type(1)=='c' && strcmp(model.latent,'no'),
    Y = 2*(Y>0)-1;
end

%
% postprocessing...
%
if model.preprocess(1)=='p' && ~(model.type(1)=='c' && strcmp(model.latent,'yes')),
    [~,Y] = postlssvm(model,[],Y);
end

%
% decode if multiclass
%
if model.type(1)=='c' && ~strcmpi(model.codetype,'none' ) && ~strcmpi(model.code,'original'),
    Y = codelssvm(model,Y);
end

%
% Simulation
%
function Y = simFct(model,X)
model.selector = ~isnan(model.ytrain);
kx = kernel_matrix(model.xtrain(model.selector, 1:model.x_dim), model.kernel_type, model.kernel_pars,X);
Y = kx'*model.alpha(model.selector,1:model.y_dim)+ones(size(kx,2),1)*model.b(:,1:model.y_dim);


function [Yt,Yl] = simmultidimoutput(model, Xt, Y,n)
%
% what to do if output multimensional?
%

Yt = []; Yl = [];
for d=1:model.y_dim,
    eval('gam = model.gam(:,d);','gam = model.gam;');
    eval('sig2 = model.kernel_pars(:,d);','sig2 = model.kernel_pars;');
    eval('kernel = model.kernel_type{d};','kernel=model.kernel_type;');
    % not yet timeseries nor NARX
    [Ytn Yln] = simlssvm({model.xtrain, model.ytrain(:,d),model.type,gam,sig2,kernel,'original'},{model.alpha(:,d),model.b(d)},Xt);
    
    Yt = [Yt Ytn]; Yl = [Yl Yln];
end

% postprocessing...
if model.preprocess(1)=='p' && ~(model.type(1)=='c' && strcmp(model.latent,'yes')),
    [~,Yt] = postlssvm(model,[],Yt);
end

% decode if multiclass
if model.type(1)=='c' && ~strcmpi(model.codetype,'none' ) && ~strcmpi(model.code,'original'),
    Yt = codelssvm(model,Yt);
end
