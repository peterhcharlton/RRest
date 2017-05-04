function model = codelssvm(model,Yt,Ytt)
% This function is only for intern LS-SVMlab use. For extern coding of the classes, use the functions 'code'.
%
% Function to decode & encode the responses of the classification model
% before applying the LS-SVM.
%
% Firstly, the appropriate functions has to be set:
%   >> model = changelssvm(model,'codetype','code_MOC');
%   >> model = changelssvm(model,'codedist_fct','codedist_bay');
% The corresponding encoding is invoked by
%   >> model = codelssvm(model)
% 
% The 2nd argument is decoded, 
%
%   Y = codelssvm(model,Y)
%
% The 3th argument is encoded, 
% 
%   Y = codelssvm(model,[],Ytt)'
%
% By default, a one dimensional categorical coding of the
% (multiclass) labels is assumed. 
% 
% ENCODE OPTIONS in the model:
%           codetype: used coding for multiclass classification;
%               code: status of the coding {'original','encoded'};
%           codetype: used coding for multiclass classification or 'none';
%       codedist_fct: function used to calculate to which class a
%                     coded result belongs;
%      codetype_args: arguments of the codetype function;
%      codedist_args: arguments of the codedist function;
%          codebook2: codebook of the new coding
%          codebook1: codebook of the original coding
%
% see also:
%      code, trainlssvm, simlssvm, code_OneVsAll, code_OneVsOne,
%      code_cat, codedist_hamming 

% (c) SCD-KULeuven, rights & help @ http://www.esat.kuleuven.ac.be/sista/lssvmlab



% default args
eval('model.codedist_args;','model.codedist_args = {};');
eval('model.codetype_args;','model.codetype_args = {};');
eval('model.codedist_fct;','model.codedist_fct = ''codedist_hamming'';');



%
% encode model
% a new rescaling is looked up in case of preprocessing
%
%
if nargin==1,
  % CLASSIFICATION %
  if model.type(1)=='c' && model.code(1)~='o', % not 'original' 
			
    if model.code(1)=='c', % 'changed'				      
            
      % check dimension
      eval('if model.y_dim~=size(model.codebook2,1), warning(''Y different dimension than original code;'');  end',...
	   'if model.y_dim~=1, warning(''Uncoded Y needs to be dimension 1;'');  end');

      % back to original coding
      eval('ys = code(model.ytrain, model.codebook1,{}, model.codebook2, model.codedist_fct,model.codedist_args{:});','ys = model.ytrain;');
      if ~strcmpi(model.codetype,'none'),
	% convert into new coding
	eval(['[ys, model.codebook2, model.codebook1,pre_yscheme] = ' ...
	      'code(ys, model.codetype, model.codetype_args, model.codebook1, model.codedist_fct, model.codedist_args{:});'],  ...
	     ['[ys, model.codebook2, model.codebook1,pre_yscheme] = ' ...
	      'code(ys, model.codetype, model.codetype_args, [],model.codedist_fct,model.codedist_args{:});']);
      end

      % postprocess - set code - preprocess
      prepro = model.preprocess; model = postlssvm(model);
      model.pre_yscheme = pre_yscheme;
      model.ytrain = ys;
      model.y_dim = size(ys,2);
      model.code = 'encoded';  
      model = changelssvm(model,'preprocess',prepro); model = prelssvm(model);
    end;
  end

  
%
% decode signal
%
elseif nargin==2,
    
  % CLASSIFICATION       %
  if model.type(1)=='c' & model.code(1)~='o', % not 'original' 
    
    eval('model.codebook1; model.codebook2;','model = trainlssvm(model);');
    if ~strcmpi(model.codetype,'none'),
      eval('model = code(Yt,model.codebook1, {}, model.codebook2, model.codedist_fct, model.codedist_args);',...
	   'model = code(Yt,model.codebook1, {}, model.codebook2);');
    else
      model=Yt;
    end
  else
    model=Yt;
  end
  
  
%
% encode signal
%
elseif nargin==3,
  
  % CLASSIFICATION       %
  if model.type(1)=='c' & model.code(1)~='o', % not 'original' 

    eval('model.codebook1; model.codebook2;','model = trainlssvm(model);');
    if ~strcmpi(model.codetype,'none'),
       eval('model = code(Ytt,model.codebook2, {}, model.codebook1, model.codedist_fct, model.codedist_args);',...
	    'model = code(Ytt,model.codebook2, {}, model.codebook1);');
     else
       model=Ytt;
     end
  else
    model=Ytt;
  end

end
  

