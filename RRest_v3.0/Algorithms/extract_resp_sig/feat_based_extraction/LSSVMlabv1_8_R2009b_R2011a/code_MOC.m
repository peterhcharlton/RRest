function [codebook,scheme] = code_MOC(m)
% Generate the codebook for multiclass classification with Minimum Output encoding.
%
% >> codebook = code_MOC(m)
%
%  see also:
%    code

% Copyright (c) 2002,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.ac.be/sista/lssvmlab

nb = ceil(log2(m));
codebook = -ones(nb,m);
for i=1:m,
  code = str2num(num2str(dec2bin(i-1)')).*2-1;
  codebook((nb-length(code)+1):nb,i) = code;
end

% output forat, where 'b' stands for binary discriminator
scheme = []; for i=1:nb, scheme = [scheme 'b']; end