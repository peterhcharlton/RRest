function [codebook,scheme] = code_OneVsOne(m)
% Generate the codebook for multiclass classification with One-Versus-One encoding.
%
% codebook = code_OneVsOne(m)
%
%  see also:
%    codelssvm

% Copyright (c) 2002,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


nb = m*(m-1)/2;
codebook = NaN*zeros(nb,m);
t=1;
for i=1:m-1,
  for j=i+1:m,
    codebook(t,i) = 1;
    codebook(t,j) = -1;
    t=t+1;
  end
end

% output format, where 'b' stands for binary discriminator
scheme = []; for i=1:nb, scheme = [scheme 'b']; end