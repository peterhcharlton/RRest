function x = latticeseq_b2(s, n, zvec)
% function x = latticeseq_b2(s, n)
%
% Generate the points from a lattice sequence in base 2 in radical inverse
% ordering.
%
% Note: the generator fist has to be initialized, e.g., by
%   latticeseq_b2 init0
% Afterwards one can obtain n s-dimensional samples in an [s x n] array as
%   x = latticeseq_b2(s, n);
%
% The default lattice sequence in this file is the ``universal'' lattice
% sequence optimized for an unanchored Sobolev space with order-2 weights from
% the paper by Cools, Kuo and Nuyens in SIAM SISC (2006).
% This sequence is ``universal'' in the sense that the weights do not really
% matter, as such, there are no good or bad weights...
% For this lattice sequence the maximum number of dimensions is 250 and the
% maximum number of points is 2^20 = 1048576. The worst-case error was
% optimized starting from 2^10 = 1024 points (this means it may perform badly
% for less than 1024 samples).
%
% Inputs:
%   s       the number of dimensions per vector
%   n       the number of samples you want
% Output:
%   x       an array of sample vectors, size [s x n]
%
% Usage:
%
%   1. Initialize the generator with the default lattice sequence
%
%     latticeseq_b2 init0
%
%   or initialize the generator with a user supplied generating vector for n
%   points
%
%     latticeseq_b2('init0', n, zvec)
%
%   Valid options for intialization are:
%     init0       the lattice sequence as is, the first point is the zero vector
%     init1       the first point is changed into an all ones vector
%     initskip    skip the first point
%
%   2. Generate the next n s-vectors of the sequence, returning an array of
%   dimensions s x n:
%
%     P = latticeseq_b2(s, n)
%
%   3. Ask the state of the current generator with
%
%     latticeseq_b2 state
%
% Example:
%   latticeseq_b2('init0');        % Initialize the generator.
%   x1 = latticeseq_b2(10, 1024);  % This gives an array of 1024 10-dimensional samples.
%   x2 = latticeseq_b2(10, 1024);  % This will give you the next 1024 10-dimensional samples.
%   latticeseq_b2('init0');        % This resets the state of the generator.
%   x3 = latticeseq_b2(10, 1024);  % This will give the same samples as in x1.
%
% (C) 2005, 2010, Dirk Nuyens, Department of Computer Science, K.U.Leuven, Belgium

persistent z smax k N b m initmode recipd maxbit

if ischar(s) && strncmp(s, 'init', 4)
    k = uint32(0);
    maxbit = 32;
    recipd = pow2(-maxbit);
    if strcmp(s, 'init0')
        k = 0;
        initmode = 0;
    elseif strcmp(s, 'init1')
        k = 0;
        initmode = 1;
    elseif strcmp(s, 'initskip')
        initmode = -1;
        k = 1;
    else
        error('I only know about ''init0'', ''init1'' and ''initskip'', what are you talking about?');
    end;
    if (nargin == 3)
        % z is given by the user
        z = zvec;
        N = n;
    else
        z = [ 1 182667 469891 498753 110745 446247 250185 118627 245333 283199 ...
              408519 391023 246327 126539 399185 461527 300343 69681 516695 436179 ...
              106383 238523 413283 70841 47719 300129 113029 123925 410745 211325 ...
              17489 511893 40767 186077 519471 255369 101819 243573 66189 152143 ...
              503455 113217 132603 463967 297717 157383 224015 502917 36237 94049 ...
              170665 79397 123963 223451 323871 303633 98567 318855 494245 477137 ...
              177975 64483 26695 88779 94497 239429 381007 110205 339157 73397 ...
              407559 181791 442675 301397 32569 147737 189949 138655 350241 63371 ...
              511925 515861 434045 383435 249187 492723 479195 84589 99703 239831 ...
              269423 182241 61063 130789 143095 471209 139019 172565 487045 304803 ...
              45669 380427 19547 425593 337729 237863 428453 291699 238587 110653 ...
              196113 465711 141583 224183 266671 169063 317617 68143 291637 263355 ...
              427191 200211 365773 254701 368663 248047 209221 279201 323179 80217 ...
              122791 316633 118515 14253 129509 410941 402601 511437 10469 366469 ...
              463959 442841 54641 44167 19703 209585 69037 33317 433373 55879 ...
              245295 10905 468881 128617 417919 45067 442243 359529 51109 290275 ...
              168691 212061 217775 405485 313395 256763 152537 326437 332981 406755 ...
              423147 412621 362019 279679 169189 107405 251851 5413 316095 247945 ...
              422489 2555 282267 121027 369319 204587 445191 337315 322505 388411 ...
              102961 506099 399801 254381 452545 309001 147013 507865 32283 320511 ...
              264647 417965 227069 341461 466581 386241 494585 201479 151243 481337 ...
              68195 75401 58359 448107 459499 9873 365117 350845 181873 7917 436695 ...
              43899 348367 423927 437399 385089 21693 268793 49257 250211 125071 ...
              341631 310163 94631 108795 21175 142847 383599 71105 65989 446433 ...
              177457 107311 295679 442763 40729 322721 420175 430359 480757 ]';
        N = pow2(20);
    end;
    smax = length(z);
    b = 2;
    m = ceil(log(N)/log(b));
    return;
elseif ischar(s) && strcmp(s, 'state')
    x.z = z;
    x.smax = smax;
    x.k = k;
    x.N = N;
    x.b = b;
    x.m = m;
    x.initmode = initmode;
    return;
end;

if ((k + n) > N) || (s > smax)
    error(sprintf('Can only generate %d lattice points in %d dimensions', N, smax));
end;

x = zeros(s, n);

if (k == 0) && (initmode == 0)
    x(:, 1) = 0; si = 2; k = k + 1;
elseif (k == 0) && (initmode == 1)
    x(:, 1) = 1; si = 2; k = k + 1;
else
   si = 1;
end;

for i=si:n
    rk = bitreverse32(k);
    x(:, i) = mod(double(rk) * recipd * z(1:s), 1);
    k = k + 1;
end;