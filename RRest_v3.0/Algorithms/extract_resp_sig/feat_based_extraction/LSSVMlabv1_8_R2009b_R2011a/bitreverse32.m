function v = bitreverse32(k)
% function v = bitreverse32(k)
%
% Reverse the bits of k.
%
% Input:
%   k       a 32 bit unsigned integer, or an array of such integers,
%           note: the input value is automatically converted to this type
% Output:
%   v       a 32 bit unsigned integer with the bits of the input in reverse
%           order, e.g., the MSB is now the LSB and vica versa
%
% See Stanford bit hacks: http://graphics.stanford.edu/~seander/bithacks.html
%
% (w) 2010, Dirk Nuyens, Department of Computer Science, KULeuven, Belgium

v = uint32(k);
% swap odd and even bits
%v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
v = bitxor( bitand( bitshift(v, -1) , 1431655765 ) , ...
            bitshift( bitand(v, 1431655765) , 1 ) );
% swap consecutive pairs
%v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
v = bitxor( bitand( bitshift(v, -2), 858993459 ) , ...
            bitshift( bitand(v, 858993459), 2 ) );
% swap nibbles ... 
%v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
v = bitxor( bitand( bitshift(v, -4), 252645135 ) , ...
            bitshift( bitand(v, 252645135), 4 ) );
% swap bytes
%v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
v = bitxor( bitand( bitshift(v, -8), 16711935 ) , ...
            bitshift( bitand(v, 16711935), 8 ) );
% swap 2-byte long pairs
%v = ( v >> 16             ) | ( v               << 16);
v = bitxor( bitshift(v, -16) , ...
            bitshift(v, 16) );