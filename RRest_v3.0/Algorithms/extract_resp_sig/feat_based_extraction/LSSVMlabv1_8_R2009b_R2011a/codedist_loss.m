function dist = codedist_loss(C1,C2,Ylat,loss)
% compute the distance using a loss function metric between 'C1' and  'C2'
%
%    dist = codedist_loss(C1,C2,Ylatent)
%    dist = codedist_loss(C1,C2,Ylatent, loss_fct)
%
% 'C1' contains the result code, 
% 'C2' contains the codebooks prototype. An infine small number 'eps'  represents the don't care
% 'loss_fct' is the loss function used in order to compute the loss
%            between 2 codewords. By default the sum of squares is used. 
%            One can do 'winner-takes-all' decoding by using the
%            loss function 'max'.
%
% An example:
%  >> Ye = [1 1; 1 1; -1 -1; 1 -1];
%  >> codebook = [1 2 3];
%  >> old_codebook = [1 1 -1; 1 -1 -1];
%  >> code(Ye, codebook, [],old_codebook,'codedist_loss',{'Ylatent','mse'})
%
% To use this distance measure in LS-SVMlab, the following
% procedure is to be followed, assume input data 'X' and multiclass
% output 'Y'
%
% % encode for training
% >> model = initlssvm(X,Y,'classification',gam,sig2,'preprocess','RBF_kernel');
% >> model = changelssvm(model,'codetype','code_OneVsOne');
% >> model = trainlssvm(model);
% 
% % decode for simulating
% >> [Yhamming, Ylatent] = simlssvm(model,Xt); 
% >> model = changelssvm(model,'codedist_fct','codedist_loss');
% >> model = changelssvm(model,'codedist_args',{Ylatent,'sse'});
% >>   Yt  = simlssvm(model,Xt);
%
% see also:
%   bay_modoutClass, codedist_hamming, code_ECOC

% (c) SCD-KULeuven, rights & help @ http://www.esat.kuleuven.ac.be/sista/lssvmlab


if nargin<3,
  warning('The latent variables should be used, proceeding with the binary classifiers...');
  Ylat = C1;
end
eval('loss;','loss=''sse'';');

[nb,nbin] = size(Ylat);
[~,dim] = size(C2);
dist = zeros(nb,dim);

for d = 1:dim,
  for n= 1:nb,
    nondontcare = find(Ylat(n,:)~=eps & C2(:,d)'~=eps);
    dist(n,d) = feval(loss, Ylat(n,nondontcare),C2(nondontcare,d)');
  end
end
dist


function l = sse(X,Y)
l = sum(sum((X-Y).^2));

function l = winnertakesall(X,Y)
p = find(Y>0);
l = max(X(P));