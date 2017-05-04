function dist = codedist_bay(C1,C2,Py)
% compute the distance using a Bayesian metric between 'C1' and  'C2'
%
%    dist = codedist_bay(C1,C2,{Py})
%
% 'C1' contains the result code, 
% 'C2' contains the codebooks prototype. '0' or  an infine small number 'eps'  represents the don't care
% 'Py' is a matrix containing the moderated outputs of the binary classifiers
%
%  
% An example:
%  >> Ye = [1 1; 1 1];
%  >> codebook = [1 2 3];
%  >> old_codebook = [1 1 -1; 1 -1 -1];
%  >> Py = [.9 -.1; -.1 .13];
%  >> code(Ye, codebook, [],old_codebook,'codedist_bay',{Py})
% in this call, 'Ye' is not used explicitly.
%
% To use this distance measure in LS-SVMlab, the following
% procedure is to be followed, assume input data 'X' and multiclass
% output 'Y':
%
%  >> [Ycode,codebook,old_codebook] = code(Y,'code_MOC');
%  >> [alpha,b] = trainlssvm({X,Yc,'c',gam,sig2});
%  >> Yhc = simlssvm({X,Yc,'c',gam,sig2},{alpha,b},Xt);
%
%  The moderated output for the LS-SVM can be computed using the bayesian inference
%  framework for the LS-SVM: 
%  >> Ymod = bay_modoutClass(model,Xt);
%  >> Yh = code(Yhc,old_codebook,[],codebook,'codedist_bay',{Ymod});
%  
% see also:
%   bay_modoutClass, codedist_hamming, code_ECOC

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

% % encode for training
% >> model = initlssvm(X,Y,'classification',gam,sig2,'preprocess','RBF_kernel');
% >> model = changelssvm(model,'codetype','code_MOC');
% >> model = changelssvm(model,'codedist_fct','codedist_hamming');
% >> model = trainlssvm(model);
% 
% % decode for simulating
% >> model = changelssvm(model,'codedist_fct','codedist_bay');
% >> model = changelssvm(model,'codedist_args',{bay_modoutClass(model,Xt)});
% >>   Yt  = simlssvm(model,Xt);


if nargin<3,
  error(['moderated output needed as function arguments' ...
	 '(for LS-SVM, model.codedist_args =' ...
	 ' bay_modoutClass(model,X)).']);
end


[nb,nbin] = size(Py);
[~,dim] = size(C2);
dist = zeros(nb,dim);

for d = 1:dim,
  for n= 1:nb,
    dist(n,d) = sum((1-Py(n,:).*C2(:,d)'))-sum(C2(:,d)==eps);
  end
end


