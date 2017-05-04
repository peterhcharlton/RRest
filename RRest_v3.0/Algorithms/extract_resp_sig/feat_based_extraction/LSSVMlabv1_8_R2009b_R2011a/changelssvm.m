function model = changelssvm(model,option, value)
% Change a field of the object oriented representation of the LS-SVM
%
%
% The different options of the fields are given in following table:
% 
%     1. General options representing the kind of model:
% 
%           type: 'classifier' ,'function estimation' 
% implementation: 'CMEX' ,'CFILE' ,'MATLAB' 
%         status: Status of this model ('trained'  or 'changed' )
%          alpha: Support values of the trained LS-SVM model
%              b: Bias term of the trained LS-SVM model
%       duration: Number of seconds the training lasts
%         latent: Returning latent variables ('no' ,'yes' ) 
%       x_delays: Number of delays of eXogeneous variables (by default 0 )
%       y_delays: Number of delays of responses (by default 0 )
%          steps: Number of steps to predict (by default 1 )
%            gam: Regularisation parameter
%    kernel_type: Kernel function
%    kernel_pars: Extra parameters of the kernel function
% 
%
%     2. Fields used to specify the used training data:
% 
%    x_dim: Dimension of input space
%    y_dim: Dimension of responses
%  nb_data: Number of training data
%   xtrain: (preprocessed) inputs of training data
%   ytrain: (preprocessed,coded) outputs of training data
% selector: Indexes of training data effectively used during training
% 
%
%     3. Options used in the Conjugate Gradient (CG) algorithm:
% 
%     cga_max_itr: Maximum number of iterations in CG
%         cga_eps: Stopcriterium of CG, largest allowed error
%    cga_fi_bound: Stopcriterium of CG, smallest allowed improvement
%        cga_show: Show the results of the CG algorithm (1 or 0)
% cga_startvalues: Starting values of the CG algorithm
% 
%
%     4. Fields with the information for pre- and post-processing (only given if appropriate):
% 
%  preprocess: 'preprocess'  or 'original' 
%     schemed: Status of the preprocessing 
%               ('coded' ,'original'  or 'schemed' )
% pre_xscheme: Scheme used for preprocessing the input data
% pre_yscheme: Scheme used for preprocessing the output data
%   pre_xmean: Mean of the input data
%    pre_xstd: Standard deviation of the input data
%   pre_ymean: Mean of the responses
%    pre_ystd: Standard deviation of the reponses
% 
%
%     5. The specifications of the used encoding (only given if appropriate):
% 
%          code: Status of the coding 
%                 ('original' ,'changed'  or 'encoded')
%      codetype: Used function for constructing the encoding 
%                  for multiclass classification (by default 'none')
% codetype_args: Arguments of the codetype function
%  codedist_fct: Function used to calculate to which class a
%                coded result belongs
% codedist_args: Arguments of the codedist function
%     codebook2: Codebook of the new coding
%     codebook1: Codebook of the original coding
% 
% Full syntax
% 
% >> model = changelssvm(model, field, value)
% 
%       Outputs    
%         model(*) : Obtained object oriented representation of the LS-SVM model
%       Inputs    
%         model    : Original object oriented representation of the LS-SVM model
%         field    : Field of the model one wants to change (e.g. 'preprocess')
%         value    : New value of the field of the model one wants to change
% 
% See also:
%   trainlssvm, initlssvm, simlssvm, plotlssvm.

% Copyright (c) 2010,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.ac.be/sista/lssvmlab


%
% alias sigma^2
%
if (strcmpi(option,'sig2')) option = 'kernel_pars'; end

%
% selector -> nb_data
% nb_data -> selector
%
if strcmp(option,'selector'),
  model.nb_data = length(value);
end
if strcmp(option,'nb_data'),
  model.selector = 1:value;
end

%
% xtrain
%
if strcmp(option,'xtrain'), 
  [nb,model.x_dim] = size(value);
  model.nb_data = nb;%min(nb,model.nb_data);
  model.selector = 1:model.nb_data;
  if length(model.gam)>model.y_dim & length(model.gam)~=size(value,1),
    warning('Discarting  different gamma''s...');
    model.gam = max(model.gam); 
  end
    eval('value=prelssvm(model,value);',...
       'warning(''new trainings inputdata not comform with used preprocessing'');');
end


%
% ytrain
%
if strcmp(option,'ytrain'),
  if size(value,2)~=size(model.ytrain,2),
    model.y_dim = size(value,2);
  end
  eval('value = codelssvm(model,[],value);',...
       'warning(''new trainings outputdata not comform with used encoding;'');');
  eval('[ff,value] = prelssvm(model,[],value);',...
       'warning(''new trainings outputdata not comform with used preprocessing;'');');
  [nb,model.y_dim] = size(value);
  model.nb_data = min(nb,model.nb_data);
  model.selector = 1:model.nb_data;
end



%
% switch between preprocessing - original data
% model.prestatus = {'changed','ok'}
%
if (strcmpi(option,'preprocess')) & model.preprocess(1)~=value(1),
  model.prestatus = 'changed'; 
end    
      
      


%
% change coding
%
if strcmpi(option,'codetype') | strcmpi(option,'codebook2') | ...
      strcmpi(option, 'codeargs') | strcmpi(option, 'codedistfct'),
  model.code = 'changed';
elseif  strcmpi(option,'codebook1'),
  warning('change original format of the classifier; the toolbox will be unable to return results in the original format');
end


%
% final change
%
eval(['old_value = model.' lower(option) ';'],'old_value=[];');
eval(['model.' lower(option) '=value;']);

if (isempty(value) | isempty(old_value)),
  different = 1;
else
  eval('different = any(old_value~=value);','different=1;');
end

if different & ~strcmpi(option,'implementation'),
  model.status = 'changed';
end

