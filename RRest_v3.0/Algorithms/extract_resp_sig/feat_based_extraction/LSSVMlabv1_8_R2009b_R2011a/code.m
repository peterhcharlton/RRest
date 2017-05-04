function [nsignals, codebook, oldcodebook, scheme] = code(signals,codetype,codetype_args,oldcodebook,fctdist,fctdist_args)
% Encode and decode a multi-class classification task into multiple binary classifiers
%
% >> Yc = code(Y, codebook)
%
% The coding is defined by the codebook. The codebook is
% represented by a matrix where the columns represent all different
% classes and the rows indicate the result of the binary
% classifiers. An example is given: the 3 classes with original
% labels [1 2 3] can be encoded in the following codebook (using Minimal Output Encoding):
%
% >> codebook
%     = [-1  -1  1;
%         1  -1  1]
% 
% For this codebook, a member of the first class is found if the
% first binary classifier is negative and the second classifier is
% positive. A don't care is represented by eps. By default it is
% assumed that the original classes are represented as different
% numerical labels. One can overrule this by passing the
% old_codebook which contains information about the old representation.
% 
% Different encoding schemes are available:
% 
%     1. Minimum Output Coding (code_MOC) 
%     2. Error Correcting Output Code (code_ECOC)
%       This coding scheme uses redundant bits. 
%     3. One versus All Coding (code_OneVsAll)
%     4. One Versus One Coding (code_OneVsOns)
% 
% Different decoding schemes are implemented:
% 
%     1. Hamming Distance (codedist_hamming) 
%     2. Bayesian Distance Measure (codedist_bay)
% 
%
% Full syntax
%
%     1. For encoding:
% 
% >> [Yc, codebook, old_codebook] = code(Y, codefct)
% >> [Yc, codebook, old_codebook] = code(Y, codefct, codefct_args)
% >> Yc = code(Y, given_codebook)
% 
%       Outputs    
%         Yc               : N x nbits encoded output classifier
%         codebook(*)      : nbits*nc matrix representing the used encoding
%         old_codebook(*)  : d*nc matrix representing the original encoding
%       Inputs    
%         Y                : N x d matrix representing the original classifier
%         codefct(*)       : Function to generate a new codebook (e.g. code_MOC)
%         codefct_args(*)  : Extra arguments for codefct
%         given_codebook(*): nbits*nc matrix representing the encoding to use
% 
%     2. For decoding:
% 
% >> Yd = code(Yc, codebook,[], old_codebook)
% >> Yd = code(Yc, codebook,[], old_codebook, codedist_fct)
% >> Yd = code(Yc, codebook,[], old_codebook, codedist_fct, codedist_args)
% 
%       Outputs    
%         Yd               : N x nc decoded output classifier
%       Inputs    
%         Y                : N x d matrix representing the original classifier
%         codebook         : d*nc matrix representing the original encoding
%         old_codebook     : bits*nc matrix representing the encoding of the given classifier
%         codedist_fct     : Function to calculate the distance between to encoded classifiers (e.g. codedist_hamming)
%         codedist_args(*) : Extra arguments of codedist_fct
% 
% 
% see also
%    code_ECOC, code_MOC, code_OneVsAll, code_OneVsOne, codedist_hamming
		    
% Copyright (c) 2010,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.ac.be/sista/lssvmlab



%
% default handling;
%
eval('fctdist(1,1);','fctdist = ''codedist_hamming'';');
eval('if isempty(oldcodebook),ss = sort(signals(:,1));  oldcodebook = ss([1;find(ss(2:end)~=ss(1:end-1))+1])'';end ',...
     'ss = sort(signals(:,1));  oldcodebook = ss([1;find(ss(2:end)~=ss(1:end-1))+1])'';');

n = size(signals,1);
mc = size(oldcodebook,2);


% codebook or codetype
% initialise the new scheme used for preprocessing
% Binary,Ctu,cAtegorical,Original
if isstr(codetype),
  eval('[codebook,scheme] = feval(codetype, mc, codetype_args{:});',...
       '[codebook,scheme] = feval(codetype, mc);');
else
  codebook = codetype;
  scheme=[]; for t=1:size(codebook,2),scheme=[scheme 'b']; end
end


%
% convert from old coding towards new coding
%
if nargin==6,
  dist = feval(fctdist, signals, oldcodebook,fctdist_args{:});
else
  dist = feval(fctdist, signals, oldcodebook);
end

for t = 1:n,
  [m,mi] = min(dist(t,:));
  m2 = min(dist(t,[1:(mi-1) (mi+1):end]));
  if m==m2, 
    nsignals(t,:) = -inf+codebook(:,mi)';
  else 
    nsignals(t,:) = codebook(:,mi)';
  end

end










