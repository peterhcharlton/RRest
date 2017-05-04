function dist = codedist_hamming(C1,C2)
% Compute the hamming distance between rows of 'C1' and  columns of 'C2'
%
% >>  distance = codedist_hamming(encoded_data, codebook);
% 
% 'encoded_data' contains the resulting codeword per row, n rows are possible
% 'codebook' contains the codebooks prototype per class as columns.
%  an infinitesimal number 'eps'  represents the don't care
%
% see also:
%   code, codelssvm, code_MOC

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

[nb,nbin] = size(C1);
[nbin,dim] = size(C2);
dist = zeros(nb,dim);
for d = 1:dim,
  for n= 1:nb,
    dist(n,d) = nbin-sum(C1(n,:)==C2(:,d)' | C1(n,:) < -10000 | C1(n,:) ==eps | C1(n,:) > 10000 | C2(:,dim)'==eps);
  end
end