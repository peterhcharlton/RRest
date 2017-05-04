function [model,Yt] = prelssvm(model,Xt,Yt)
% Preprocessing of the LS-SVM
%
% These functions should only be called by trainlssvm or by
% simlssvm. At first the preprocessing assigns a label to each in-
% and output component (c for continuous, a for categorical or b
% for binary variables). According to this label each dimension is rescaled:
% 
%     * continuous: zero mean and unit variance
%     * categorical: no preprocessing
%     * binary: labels -1 and +1
% 
% Full syntax (only using the object oriented interface):
% 
% >> model   = prelssvm(model)
% >> Xp = prelssvm(model, Xt)
% >> [empty, Yp] = prelssvm(model, [], Yt)
% >> [Xp, Yp] = prelssvm(model, Xt, Yt)
% 
%       Outputs    
%         model : Preprocessed object oriented representation of the LS-SVM model
%         Xp    : Nt x d matrix with the preprocessed inputs of the test data
%         Yp    : Nt x d matrix with the preprocessed outputs of the test data
%       Inputs    
%         model : Object oriented representation of the LS-SVM model
%         Xt    : Nt x d matrix with the inputs of the test data to preprocess
%         Yt    : Nt x d matrix with the outputs of the test data to preprocess
% 
% 
% See also:
%   postlssvm, trainlssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

if model.preprocess(1)~='p', % no 'preprocessing
  if nargin>=2, model = Xt;  end 
  return
end


% 
% what to do
% 
if model.preprocess(1)=='p', 
  eval('if model.prestatus(1)==''c'',model.prestatus=''unschemed'';end','model.prestatus=''unschemed'';');
end  


if nargin==1, % only model rescaling		
  %
  % if UNSCHEMED, redefine a rescaling
  %
  if model.prestatus(1)=='u',% 'unschemed'
    ffx =[];
    
    
    for i=1:model.x_dim,
      eval('ffx = [ffx model.pre_xscheme(i)];',...
	   'ffx = [ffx signal_type(model.xtrain(:,i),inf)];');
    end
    model.pre_xscheme = ffx;
   
    ff = [];
    for i=1:model.y_dim,
      eval('ff = [ff model.pre_yscheme(i)];',...
	   'ff = [ff signal_type(model.ytrain(:,i),model.type)];');
    end
    model.pre_yscheme = ff;
    model.prestatus='schemed';
  end
  
  %
  % execute rescaling as defined if not yet CODED
  %
  if model.prestatus(1)=='s',% 'schemed'  
    model=premodel(model); 
    model.prestatus = 'ok';
  end
  
  %
  % rescaling of the to simulate inputs
  %
elseif model.preprocess(1)=='p'
  if model.prestatus(1)=='o',%'ok' 
    eval('Yt;','Yt=[];');
    [model,Yt] = premodel(model,Xt,Yt);
  else 
    warning('model rescaling inconsistent..redo ''model=prelssvm(model);''..');
  end
end





function [type,ss] = signal_type(signal,type)
%
% determine the type of the signal,
% binary classifier ('b'), categorical classifier ('a'), or continuous
% signal ('c')
%
%
ss = sort(signal);
dif = sum(ss(2:end)~=ss(1:end-1))+1;
% binary
if dif==2,
  type = 'b';

% categorical
elseif dif<sqrt(length(signal)) || type(1)== 'c',
  type='a';

% continu
else
  type ='c';
end
  




%
% effective rescaling
%
function [model,Yt] = premodel(model,Xt,Yt)
%
%
%

if nargin==1,

  for i=1:model.x_dim,
    % CONTINUOUS VARIABLE: 
    if model.pre_xscheme(i)=='c',
      model.pre_xmean(i)=mean(model.xtrain(:,i));
      model.pre_xstd(i) = std(model.xtrain(:,i));
      model.xtrain(:,i) = pre_zmuv(model.xtrain(:,i),model.pre_xmean(i),model.pre_xstd(i));
      % CATEGORICAL VARIBALE: 
    elseif model.pre_xscheme(i)=='a',
      model.pre_xmean(i)= 0;
      model.pre_xstd(i) = 0;
      model.xtrain(:,i) = pre_cat(model.xtrain(:,i),model.pre_xmean(i),model.pre_xstd(i));
      % BINARY VARIBALE: 
    elseif model.pre_xscheme(i)=='b',      
      model.pre_xmean(i) = min(model.xtrain(:,i));
      model.pre_xstd(i) = max(model.xtrain(:,i));
      model.xtrain(:,i) = pre_bin(model.xtrain(:,i),model.pre_xmean(i),model.pre_xstd(i));
    end  
  end
  
  for i=1:model.y_dim,
    % CONTINUOUS VARIABLE: 
    if model.pre_yscheme(i)=='c',
      model.pre_ymean(i)=mean(model.ytrain(:,i),1);
      model.pre_ystd(i) = std(model.ytrain(:,i),1);
      model.ytrain(:,i) = pre_zmuv(model.ytrain(:,i),model.pre_ymean(i),model.pre_ystd(i));
    % CATEGORICAL VARIBALE: 
    elseif model.pre_yscheme(i)=='a',      
      model.pre_ymean(i)=0;
      model.pre_ystd(i) =0;
      model.ytrain(:,i) = pre_cat(model.ytrain(:,i),model.pre_ymean(i),model.pre_ystd(i));
    % BINARY VARIBALE: 
    elseif model.pre_yscheme(i)=='b',      
      model.pre_ymean(i) = min(model.ytrain(:,i));
      model.pre_ystd(i) = max(model.ytrain(:,i));
      model.ytrain(:,i) = pre_bin(model.ytrain(:,i),model.pre_ymean(i),model.pre_ystd(i));
    end  
  end

else %if nargin>1, % testdata Xt, 
  if ~isempty(Xt),
    if size(Xt,2)~=model.x_dim, warning('dimensions of Xt not compatible with dimensions of support vectors...');end
    for i=1:model.x_dim,
      % CONTINUOUS VARIABLE: 
      if model.pre_xscheme(i)=='c',
	Xt(:,i) = pre_zmuv(Xt(:,i),model.pre_xmean(i),model.pre_xstd(i));
      % CATEGORICAL VARIBALE: 
      elseif model.pre_xscheme(i)=='a',
	Xt(:,i) = pre_cat(Xt(:,i),model.pre_xmean(i),model.pre_xstd(i));
      % BINARY VARIBALE: 
      elseif model.pre_xscheme(i)=='b',      
	Xt(:,i) = pre_bin(Xt(:,i),model.pre_xmean(i),model.pre_xstd(i));
      end  
    end
  end
  
  if nargin>2 & ~isempty(Yt),
    if size(Yt,2)~=model.y_dim, warning('dimensions of Yt not compatible with dimensions of training output...');end
    for i=1:model.y_dim,
      % CONTINUOUS VARIABLE: 
      if model.pre_yscheme(i)=='c',
	Yt(:,i) = pre_zmuv(Yt(:,i),model.pre_ymean(i), model.pre_ystd(i));
      % CATEGORICAL VARIBALE: 
      elseif model.pre_yscheme(i)=='a',      
	Yt(:,i) = pre_cat(Yt(:,i),model.pre_ymean(i),model.pre_ystd(i));
      % BINARY VARIBALE: 
      elseif model.pre_yscheme(i)=='b',      
	Yt(:,i) = pre_bin(Yt(:,i),model.pre_ymean(i),model.pre_ystd(i));
      end
    end
  end
  
  % assign output
  model=Xt;
end


function X = pre_zmuv(X,mean,var)
%
% preprocessing a continuous signal; rescaling to zero mean and unit
% variance 
% 'c'
%
X = (X-mean)./var;


function X = pre_cat(X,mean,range)
%
% preprocessing a categorical signal;
% 'a'
%
X=X;


function X = pre_bin(X,min,max)
%
% preprocessing a binary signal;
% 'b'
%
if ~sum(isnan(X)) >= 1 %--> OneVsOne encoding
    n = (X==min);
    p = not(n);
    X=-1.*(n)+p;
end



