function [codebook,scheme] = code_OneVsAll(m)
% Generate the codebook for multiclass classification with One-Versus-All encoding.
%
% codebook = code_OneVsAll(m)
%
%  see also:
%    code, codedist_hamming

% (c) SCD-KULeuven, rights & help @ http://www.esat.kuleuven.be/sista/lssvmlab


codebook = eye(m).*2-1;
scheme = []; for i=1:m, scheme = [scheme 'b']; end