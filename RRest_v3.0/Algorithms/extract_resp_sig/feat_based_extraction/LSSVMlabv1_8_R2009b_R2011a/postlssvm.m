function [model,Yt] = postlssvm(model,Xt,Yt)
% Postprocessing of the LS-SVM
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
% >> model   = postssvm(model)
% >> Xp = postlssvm(model, Xt)
% >> [empty, Yp] = postlssvm(model, [], Yt)
% >> [Xp, Yp] = postlssvm(model, Xt, Yt)
% 
%       Outputs    
%         model : Preprocessed object oriented representation of the LS-SVM model
%         Xt    : Nt x d matrix with the inputs of the test data to preprocess
%         Yt    : Nt x d matrix with the outputs of the test data to preprocess
%       Inputs    
%         model : Object oriented representation of the LS-SVM model
%         Xp    : Nt x d matrix with the preprocessed inputs of the test data
%         Yp    : Nt x d matrix with the preprocessed outputs of the test data
% 
% 
% See also:
%   prelssvm, trainlssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

%
% test of postprocessing needed, and if model is properly coded before
% decoding 
%
if model.preprocess(1)~='p',
  if nargin>=2,
    % no postprocessing nor coding needed
    model=Xt; 
  end
  return
end


%
% postprocess the LS-SVM
%
if nargin==1,  
     
    %
    % execute rescaling as defined 
    %
    if (model.prestatus(1)=='o' & model.preprocess(1)=='p') | ... % 'preprocess' &'ok'
       (model.prestatus(1)=='c' & model.preprocess(1)=='o'),      % 'original' &'changed'
      model=postmodel(model);   
      model.preprocess = 'original'; 
    end
    model.prestatus='ok';
    
   

    %
    % rescaling of the to simulate inputs
    %
else
  eval('Yt;','Yt=[];');
  [model,Yt] = postmodel(model,Xt,Yt);
end





function [model,Yt] = postmodel(model,Xt,Yt)
%
% ' [Xt,Yt] = postmodel(model,Xt,Yt)'
% ' [model] = postmodel(model)'
%

if nargin==1,
  
  for i=1:model.x_dim,
    % CONTINU VARIABLE: 
    if model.pre_xscheme(i)=='c',
      model.xtrain(:,i) = post_zmuv(model.xtrain(:,i),model.pre_xmean(i),model.pre_xstd(i));
    % CATHEGORICAL VARIBALE: 
    elseif model.pre_xscheme(i)=='a',
      model.xtrain(:,i) = post_cat(model.xtrain(:,i),model.pre_xmean(i),model.pre_xstd(i));
    % BINARY VARIBALE: 
    elseif model.pre_xscheme(i)=='b',
      model.xtrain(:,i) = post_bin(model.xtrain(:,i),model.pre_xmean(i),model.pre_xstd(i));
    end
  end
  
  for i=1:model.y_dim,
    % CONTINU VARIABLE: 
    if model.pre_yscheme(i)=='c',
      model.ytrain(:,i) = post_zmuv(model.ytrain(:,i),model.pre_ymean(i),model.pre_ystd(i));
    % CATHEGORICAL VARIBALE: 
    elseif model.pre_yscheme(i)=='a',      
      model.ytrain(:,i) = post_cat(model.ytrain(:,i),model.pre_ymean(i),model.pre_ystd(i));
    % BINARY VARIBALE: 
    elseif model.pre_yscheme(i)=='b',
      model.ytrain(:,i) = post_bin(model.ytrain(:,i),model.pre_ymean(i),model.pre_ystd(i));
    end  
  end

else
  
  if nargin>1, % testdayta Xt,
    if ~isempty(Xt),
      if size(Xt,2)~=model.x_dim, warning('dimensions of Xt not compatible with dimensions of supprt vectors...');end
      for i=1:model.x_dim,
	% CONTINU VARIABLE: 
	if model.pre_xscheme(i)=='c',
	  Xt(:,i) = post_zmuv(Xt(:,i),model.pre_xmean(i),model.pre_xstd(i));
	  % CATHEGORICAL VARIBALE: 
	elseif model.pre_xscheme(i)=='a',
	  Xt(:,i) = post_cat(Xt(:,i),model.pre_xmean(i),model.pre_xstd(i));
	% BINARY VARIBALE: 
	elseif model.pre_xscheme(i)=='b',
	  Xt(:,i) = post_bin(Xt(:,i),model.pre_xmean(i),model.pre_xstd(i));
	end
      end
    end
    
    if nargin>2 & ~isempty(Yt),
      if size(Yt,2)~=model.y_dim, warning('dimensions of Yt not compatible with dimensions of supprt vectors...');end
      for i=1:model.y_dim,
	% CONTINU VARIABLE: 
	if model.pre_yscheme(i)=='c',
	  Yt(:,i) = post_zmuv(Yt(:,i),model.pre_ymean(i), model.pre_ystd(i));
	  % CATHEGORICAL VARIBALE: 
	elseif model.pre_yscheme(i)=='a',      
	  Yt(:,i) = post_cat(Yt(:,i),model.pre_ymean(i),model.pre_ystd(i));
	% BINARY VARIBALE: 
	elseif model.pre_yscheme(i)=='b',
	  Yt(:,i) = post_bin(Yt(:,i),model.pre_ymean(i),model.pre_ystd(i));
	end
      end
    end
    model = Xt;
  end
end



function X = post_zmuv(X,mean,var)
%
% postprocessing a continu signal; rescaling to zero mean and unit
% variance 
% 'c'
%
X = X.*var+mean;


function X = post_cat(X,mean,range)
%
% postprocessing a cathegorical signal, rescaling to -1:1;
% 'a'
%
X = X;


function X = post_bin(X,min,max)
%
% postprocessing a binary signal, rescaling to -1:1;
% 'a'
%
X = min.*(X<=0)+max.*(X>0);
